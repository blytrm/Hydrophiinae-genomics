//! Parallel cyclic null for V2R cluster × LTR coverage.
//!
//! Matches `regen_empirical_nulls.R`:
//!   coverage = Σ|LTR ∩ cluster| / Σ|cluster|
//!   null     = circularly rotate phase (here: shift clusters by −offset)
//!
//! Speed tricks:
//!   1. Shift the *small* cluster set, not the huge LTR track (same geometry).
//!   2. Pre-reduce / sort LTR once; binary-search overlaps.
//!   3. Rayon parallelises the permutations across CPU cores.

use anyhow::{bail, Context, Result};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// One genomic interval, 1-based inclusive (GenomicRanges style).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Interval {
    start: i64,
    end: i64,
}

impl Interval {
    #[inline]
    fn width(self) -> i64 {
        self.end - self.start + 1
    }
}

/// Merge overlapping / abutting intervals (GenomicRanges::reduce).
fn reduce_intervals(mut ivs: Vec<Interval>) -> Vec<Interval> {
    if ivs.is_empty() {
        return ivs;
    }
    ivs.sort_unstable_by_key(|x| x.start);
    let mut out = Vec::with_capacity(ivs.len());
    let mut cur = ivs[0];
    for nxt in ivs.into_iter().skip(1) {
        if nxt.start <= cur.end + 1 {
            cur.end = cur.end.max(nxt.end);
        } else {
            out.push(cur);
            cur = nxt;
        }
    }
    out.push(cur);
    out
}

/// Merge if gap < `min_gap` (GenomicRanges::reduce(min.gapwidth = G)).
fn reduce_with_gap(mut ivs: Vec<Interval>, min_gap: i64) -> Vec<Interval> {
    if ivs.is_empty() {
        return ivs;
    }
    ivs.sort_unstable_by_key(|x| x.start);
    let mut out = Vec::with_capacity(ivs.len());
    let mut cur = ivs[0];
    for nxt in ivs.into_iter().skip(1) {
        if nxt.start <= cur.end + min_gap {
            cur.end = cur.end.max(nxt.end);
        } else {
            out.push(cur);
            cur = nxt;
        }
    }
    out.push(cur);
    out
}

/// Circularly rotate intervals on a chromosome of length `chr_len` by `offset`.
/// Identical to the R `cyclic_shift` in regen_empirical_nulls.R.
fn cyclic_shift(ivs: &[Interval], chr_len: i64, offset: i64) -> Vec<Interval> {
    if ivs.is_empty() {
        return Vec::new();
    }
    let mut parts = Vec::with_capacity(ivs.len() * 2);
    for iv in ivs {
        let w = iv.width();
        let ns = ((iv.start + offset - 1).rem_euclid(chr_len)) + 1;
        let ne = ns + w - 1;
        if ne <= chr_len {
            parts.push(Interval {
                start: ns,
                end: ne,
            });
        } else {
            // wraps past the end → two pieces
            parts.push(Interval {
                start: ns,
                end: chr_len,
            });
            parts.push(Interval {
                start: 1,
                end: ne - chr_len,
            });
        }
    }
    reduce_intervals(parts)
}

/// Σ |A ∩ B| / Σ |A|.  A is reduced; B must already be sorted non-overlapping.
fn coverage_fraction(a: &[Interval], b: &[Interval]) -> f64 {
    let a = reduce_intervals(a.to_vec());
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let abp: i64 = a.iter().map(|x| x.width()).sum();
    if abp == 0 {
        return 0.0;
    }

    // b starts / ends for binary search
    let b_starts: Vec<i64> = b.iter().map(|x| x.start).collect();
    let b_ends: Vec<i64> = b.iter().map(|x| x.end).collect();

    let mut covered: i64 = 0;
    for iv in &a {
        // last index with b.start <= iv.end
        let i0 = b_starts.partition_point(|&s| s <= iv.end);
        // first index with b.end >= iv.start
        let j0 = b_ends.partition_point(|&e| e < iv.start);
        if j0 >= i0 {
            continue;
        }
        for k in j0..i0 {
            let ov_s = iv.start.max(b[k].start);
            let ov_e = iv.end.min(b[k].end);
            let w = ov_e - ov_s + 1;
            if w > 0 {
                covered += w;
            }
        }
    }
    covered as f64 / abp as f64
}

/// One chromosome's cyclic null (parallel over permutations).
fn cyclic_null_parallel(
    clusters: &[Interval],
    ltr: &[Interval],
    chr_len: i64,
    ntimes: usize,
    seed: u64,
) -> (f64, Vec<f64>) {
    let observed = coverage_fraction(clusters, ltr);

    // Draw offsets with a serial RNG first (reproducible), then eval in parallel.
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let offsets: Vec<i64> = (0..ntimes)
        .map(|_| rng.gen_range(1..=chr_len))
        .collect();

    let values: Vec<f64> = offsets
        .par_iter()
        .map(|&off| {
            // cov(A, shift(B, off)) == cov(shift(A, −off), B)
            let neg = (-off).rem_euclid(chr_len);
            let shifted = cyclic_shift(clusters, chr_len, neg);
            coverage_fraction(&shifted, ltr)
        })
        .collect();

    (observed, values)
}

fn summarise(obs: f64, vals: &[f64]) -> (f64, f64, f64, f64, f64) {
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    let sd = var.sqrt();
    let fold = obs / mean;
    let z = (obs - mean) / sd;
    let ge = vals.iter().filter(|&&v| v >= obs).count();
    let p = (ge as f64 + 1.0) / (n + 1.0);
    (mean, sd, fold, z, p)
}

// ───────────────────────── I/O helpers ─────────────────────────

fn read_fai(path: &Path) -> Result<HashMap<String, i64>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut map = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let chr = parts.next().context("fai chr")?.to_string();
        let len: i64 = parts.next().context("fai len")?.parse()?;
        map.insert(chr, len);
    }
    Ok(map)
}

/// BED0 → 1-based inclusive intervals, optionally restricted to `keep`.
fn read_bed0(path: &Path, keep: &HashSet<String>) -> Result<HashMap<String, Vec<Interval>>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out: HashMap<String, Vec<Interval>> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let mut chr = parts.next().context("bed chr")?.to_string();
        if chr == "ch6" {
            chr = "chZ".into();
        }
        if !keep.contains(&chr) {
            continue;
        }
        let start0: i64 = parts.next().context("bed start")?.parse()?;
        let end: i64 = parts.next().context("bed end")?.parse()?;
        out.entry(chr).or_default().push(Interval {
            start: start0 + 1,
            end,
        });
    }
    Ok(out)
}

fn intervals_overlap(a: Interval, b_sorted: &[Interval]) -> bool {
    let starts: Vec<i64> = b_sorted.iter().map(|x| x.start).collect();
    let ends: Vec<i64> = b_sorted.iter().map(|x| x.end).collect();
    let i0 = starts.partition_point(|&s| s <= a.end);
    let j0 = ends.partition_point(|&e| e < a.start);
    if j0 >= i0 {
        return false;
    }
    for k in j0..i0 {
        let ov_s = a.start.max(b_sorted[k].start);
        let ov_e = a.end.min(b_sorted[k].end);
        if ov_e >= ov_s {
            return true;
        }
    }
    false
}

fn load_te_class(path: &Path) -> Result<HashMap<String, String>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut map = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let id = parts.next().context("te id")?.to_string();
        let class = parts.next().context("te class")?.to_string();
        map.entry(id).or_insert(class); // first wins (= unique by te_id)
    }
    Ok(map)
}

fn load_ltr(
    repeats_bed: &Path,
    te_tsv: &Path,
    keep: &HashSet<String>,
) -> Result<HashMap<String, Vec<Interval>>> {
    let te = load_te_class(te_tsv)?;
    let f = File::open(repeats_bed).with_context(|| format!("open {}", repeats_bed.display()))?;
    let mut raw: HashMap<String, Vec<Interval>> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let mut chr = parts.next().context("rep chr")?.to_string();
        if chr == "ch6" {
            chr = "chZ".into();
        }
        if !keep.contains(&chr) {
            continue;
        }
        let start0: i64 = parts.next().context("rep start")?.parse()?;
        let end: i64 = parts.next().context("rep end")?.parse()?;
        let name = parts.next().context("rep name")?;
        match te.get(name) {
            Some(class) if class.starts_with("LTR/") => {
                raw.entry(chr).or_default().push(Interval {
                    start: start0 + 1,
                    end,
                });
            }
            _ => {}
        }
    }
    let mut out = HashMap::new();
    for (chr, ivs) in raw {
        out.insert(chr, reduce_intervals(ivs));
    }
    Ok(out)
}

fn default_datadir() -> PathBuf {
    PathBuf::from("/hpcfs/users/a1864358/sanders_lab/repeatdensity-atv2rs/inv/tech-rep")
}

fn default_meta() -> PathBuf {
    PathBuf::from(
        "/hpcfs/users/a1864358/sanders_lab/repeatdensity-atv2rs/final-bt/data/results/empirical_null_meta_ltr_coverage.tsv",
    )
}

/// Load R-exported offsets: columns i, offset, value
fn read_offset_table(path: &Path) -> Result<(Vec<i64>, Vec<f64>)> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut offsets = Vec::new();
    let mut values = Vec::new();
    for (li, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        if li == 0 && line.starts_with("i\t") {
            continue;
        }
        let mut parts = line.split('\t');
        let _i = parts.next().context("i")?;
        let off: i64 = parts.next().context("offset")?.parse()?;
        let val: f64 = parts.next().context("value")?.parse()?;
        offsets.push(off);
        values.push(val);
    }
    Ok((offsets, values))
}

fn load_analysis_inputs(
    datadir: &Path,
    merge_gap: i64,
) -> Result<(
    HashMap<String, i64>,
    HashMap<String, Vec<Interval>>,
    HashMap<String, Vec<Interval>>,
)> {
    let chrs = ["ch2", "chZ"];
    let keep: HashSet<String> = chrs.iter().map(|s| (*s).to_string()).collect();
    let chrlen = read_fai(&datadir.join("horn.fa.fai"))?;
    for ch in &chrs {
        if !chrlen.contains_key(*ch) {
            bail!("missing {ch} in fai");
        }
    }

    let phmm = read_bed0(&datadir.join("phmm_intact_loci_c2cZ.bed"), &keep)?;
    let tblastn = read_bed0(&datadir.join("tblastn_loci_c2cZ.bed"), &keep)?;
    let mut clusters: HashMap<String, Vec<Interval>> = HashMap::new();
    let mut n_anchor = 0usize;
    for ch in &chrs {
        let tb = reduce_intervals(tblastn.get(*ch).cloned().unwrap_or_default());
        let mut hits = Vec::new();
        if let Some(ivs) = phmm.get(*ch) {
            for &iv in ivs {
                if intervals_overlap(iv, &tb) {
                    hits.push(iv);
                }
            }
        }
        n_anchor += hits.len();
        // GenomicRanges::reduce(min.gapwidth = merge_gap)
        clusters.insert((*ch).to_string(), reduce_with_gap(hits, merge_gap));
    }
    eprintln!("anchor loci={n_anchor}  merge_gap={merge_gap}");

    let t_ltr = Instant::now();
    let ltr = load_ltr(
        &datadir.join("repeats_c2cZ.bed"),
        &datadir.join("te_classification.tsv"),
        &keep,
    )?;
    eprintln!("LTR load+reduce in {:.2?}", t_ltr.elapsed());
    Ok((chrlen, clusters, ltr))
}

/// Validate against R: same offsets → identical coverage values; observed vs notebook meta.
fn run_validate(datadir: &Path, offset_dir: &Path, meta_path: &Path) -> Result<()> {
    let (chrlen, clusters, ltr) = load_analysis_inputs(datadir, 39_810)?;
    let meta = {
        let f = File::open(meta_path).with_context(|| format!("open {}", meta_path.display()))?;
        let mut map = HashMap::new();
        for (li, line) in BufReader::new(f).lines().enumerate() {
            let line = line?;
            if li == 0 {
                continue;
            }
            let p: Vec<&str> = line.split('\t').collect();
            if p.len() < 8 || p[0] != "Cyclic" {
                continue;
            }
            // null, chr, observed, null_mean, null_sd, fold, z, p_greater
            map.insert(
                p[1].to_string(),
                (
                    p[2].parse::<f64>()?,
                    p[3].parse::<f64>()?,
                    p[7].parse::<f64>()?,
                ),
            );
        }
        map
    };

    let mut all_ok = true;
    println!(
        "chr\tn_cl\tobs_match\tperms_exact\tmax_abs_delta\tpy_p\tR_file_p\tR_meta_p\tnull_mean\tR_meta_mean\tsecs"
    );

    for ch in ["ch2", "chZ"] {
        let path = offset_dir.join(format!("cyclic_offsets_{ch}.tsv"));
        let (offsets, r_vals) = read_offset_table(&path)?;
        let cl = clusters.get(ch).context("clusters")?;
        let lt = ltr.get(ch).context("ltr")?;
        let chr_len = *chrlen.get(ch).unwrap();
        let obs = coverage_fraction(cl, lt);

        let t0 = Instant::now();
        let rust_vals: Vec<f64> = offsets
            .par_iter()
            .map(|&off| {
                let neg = (-off).rem_euclid(chr_len);
                coverage_fraction(&cyclic_shift(cl, chr_len, neg), lt)
            })
            .collect();
        let secs = t0.elapsed().as_secs_f64();

        let n = rust_vals.len();
        let mut n_exact = 0usize;
        let mut max_abs = 0.0f64;
        for (&a, &b) in rust_vals.iter().zip(r_vals.iter()) {
            let d = (a - b).abs();
            if d > max_abs {
                max_abs = d;
            }
            // match Python validation tolerance
            if d <= 1e-12 {
                n_exact += 1;
            }
        }
        let (_, _, _, _, p_rust) = summarise(obs, &rust_vals);
        let p_rfile = {
            let ge = r_vals.iter().filter(|&&v| v >= obs).count();
            (ge as f64 + 1.0) / (n as f64 + 1.0)
        };
        let mean_rust = rust_vals.iter().sum::<f64>() / n as f64;
        let (r_obs, r_mean, r_p) = meta.get(ch).copied().unwrap_or((f64::NAN, f64::NAN, f64::NAN));
        let obs_match = (obs - r_obs).abs() <= 1e-12;

        println!(
            "{ch}\t{}\t{obs_match}\t{n_exact}/{n}\t{max_abs:.3e}\t{p_rust:.6}\t{p_rfile:.6}\t{r_p:.6}\t{mean_rust:.6}\t{r_mean:.6}\t{secs:.4}",
            cl.len()
        );

        if !obs_match || n_exact != n {
            all_ok = false;
        }
        eprintln!(
            "  {ch}: obs_rust={obs:.15} obs_meta={r_obs:.15} Δ={:.3e}",
            obs - r_obs
        );
    }

    if all_ok {
        eprintln!("VALIDATION PASS: Rust == R on all matched offsets; observed matches notebook meta");
        Ok(())
    } else {
        bail!("VALIDATION FAIL: see table above");
    }
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    // Usage:
    //   cyclic_null [datadir] [ntimes] [seed] [merge_gap]
    //   cyclic_null --validate [offset_dir] [datadir] [meta.tsv]
    if args.get(1).map(|s| s.as_str()) == Some("--validate") {
        let offset_dir = PathBuf::from(args.get(2).map(|s| s.as_str()).unwrap_or("/tmp"));
        let datadir = args
            .get(3)
            .map(PathBuf::from)
            .unwrap_or_else(default_datadir);
        let meta = args.get(4).map(PathBuf::from).unwrap_or_else(default_meta);
        eprintln!(
            "validate: offsets={} datadir={} meta={}",
            offset_dir.display(),
            datadir.display(),
            meta.display()
        );
        return run_validate(&datadir, &offset_dir, &meta);
    }

    let datadir = if args.len() > 1 {
        PathBuf::from(&args[1])
    } else {
        default_datadir()
    };
    let ntimes: usize = args
        .get(2)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(2000);
    let seed: u64 = args
        .get(3)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(13);
    // default 39810 bp (~39.81 kb) — between modes of V2R nearest-gap distribution
    let merge_gap: i64 = args
        .get(4)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(39_810);

    let chrs = ["ch2", "chZ"];
    eprintln!(
        "datadir={}  ntimes={}  seed={}  merge_gap={}",
        datadir.display(),
        ntimes,
        seed,
        merge_gap
    );
    let t_all = Instant::now();
    let (chrlen, clusters, ltr) = load_analysis_inputs(&datadir, merge_gap)?;

    println!("chr\tn_clusters\tobserved\tnull_mean\tnull_sd\tfold\tz\tp_greater\tntimes\tsecs");

    for ch in &chrs {
        let cl = clusters.get(*ch).context("clusters")?;
        let lt = ltr.get(*ch).context("ltr")?;
        let chr_len = *chrlen.get(*ch).unwrap();
        eprintln!(
            "{ch}: {} clusters, {} LTR intervals, chr_len={chr_len}",
            cl.len(),
            lt.len()
        );

        let t0 = Instant::now();
        let chrom_seed = seed
            ^ (ch
                .as_bytes()
                .iter()
                .fold(0u64, |a, &b| a.wrapping_mul(31).wrapping_add(b as u64)));
        let (obs, vals) = cyclic_null_parallel(cl, lt, chr_len, ntimes, chrom_seed);
        let (mean, sd, fold, z, p) = summarise(obs, &vals);
        let secs = t0.elapsed().as_secs_f64();

        // null draws for violin plots
        let draws_path = datadir.join(format!("cyclic_draws_{ch}.tsv"));
        {
            let mut f = File::create(&draws_path)
                .with_context(|| format!("create {}", draws_path.display()))?;
            writeln!(f, "value")?;
            for v in &vals {
                writeln!(f, "{v:.15}")?;
            }
        }

        println!(
            "{ch}\t{}\t{obs:.15}\t{mean:.15}\t{sd:.15}\t{fold:.15}\t{z:.15}\t{p:.15}\t{ntimes}\t{secs:.4}",
            cl.len()
        );
    }

    eprintln!("total {:.2?}", t_all.elapsed());
    Ok(())
}
