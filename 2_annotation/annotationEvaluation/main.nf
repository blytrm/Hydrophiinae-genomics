#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// params 
params.proteins      = null          // proteins (.faa) 
params.cds           = null          //  CDS  (.fna)
params.gff3          = null          // predicted annotation
params.genome_fasta  = null          // genome (for premature-stop flag + CpG plot)
params.gff3_ref      = null          // optional reference annotation
params.busco_lineage = 'squamata'
params.compleasm_lib = null          // local compleasm/BUSCO lineage folder (-L)
params.omark_db      = null          // omamer .h5 database
params.ete_db        = null          // prebuilt ete3 NCBI taxa.sqlite (offline OMArk)
params.outdir        = 'results'

// hard requirements
[ proteins: params.proteins, cds: params.cds, gff3: params.gff3,
  genome_fasta: params.genome_fasta, omark_db: params.omark_db,
  ete_db: params.ete_db,
  compleasm_lib: params.compleasm_lib ].each { k, v ->
    if (!v) error "missing required param: --${k}"
}

def modules = "${projectDir}/modules"

// workflow
workflow {
    ch_proteins = Channel.fromPath(params.proteins, checkIfExists: true)
    ch_cds      = Channel.fromPath(params.cds,      checkIfExists: true)
    ch_gff3     = Channel.fromPath(params.gff3,     checkIfExists: true)
    ch_genome   = Channel.fromPath(params.genome_fasta, checkIfExists: true)
    ch_db       = Channel.fromPath(params.omark_db, checkIfExists: true)
    ch_ete      = Channel.fromPath(params.ete_db,   checkIfExists: true)
    ch_ref = params.gff3_ref ? Channel.fromPath(params.gff3_ref, checkIfExists: true) : Channel.empty()

    // evaluation
    PSAURON(ch_cds, ch_proteins, file("${modules}/extract_cds.py"))
    COMPLEASM(ch_proteins, params.busco_lineage, file(params.compleasm_lib),
              file("${modules}/compleasm_offline.py"))
    OMARK(ch_proteins, ch_db, ch_ete, file("${modules}/make_isoform_file.py"))
    AGAT_STATS(ch_gff3, ch_genome)
    GFFUTILS(ch_gff3, file("${modules}/gffutils_script.py"))

    // plots
    PLOT_GENELENGTH(ch_gff3,   file("${modules}/plot_genelength.py"))
    PLOT_GENEDENSITY(ch_gff3,  file("${modules}/plot_genedensity.py"))
    PLOT_CPG(ch_genome,        file("${modules}/plot_cpg.R"))

    // reference/original comparison (if --gff3_ref given)
    COMPARE_REF(ch_ref, ch_gff3, file("${modules}/compare_toRef.py"))

    // gather
    ch_reports = PSAURON.out.summary
        .mix( COMPLEASM.out.summary,
              OMARK.out.summary,
              AGAT_STATS.out.stats,
              GFFUTILS.out.txt,
              params.gff3_ref ? COMPARE_REF.out.summary : Channel.empty() )
        .collect()

    SUMMARY_REPORT(ch_reports)
}

// processes
process PSAURON {
    tag   { cds.simpleName }
    label 'psauron'
    conda 'bioconda::psauron'
    publishDir "${params.outdir}/psauron", mode: 'copy'

    input:
    path cds
    path proteins
    path extract_cds

    output:
    path "psauron_cds_${cds.simpleName}.csv"
    path "psauron_protein_${cds.simpleName}.csv"
    path "psauron_summary_${cds.simpleName}.txt", emit: summary

    script:
    def base = cds.simpleName
    """
    # Headline score: spliced CDS scored by the full (in-frame vs out-of-frame)
    # model. Transcripts carry UTRs, so cut the CDS span out first.
    python3 ${extract_cds} ${cds} cds_${base}.fna
    psauron -i cds_${base}.fna -c -o psauron_cds_${base}.csv | tee psauron_cds_${base}.log

    # Protein-mode score on the annotated amino-acid sequences (-p).
    psauron -i ${proteins} -p -c -o psauron_protein_${base}.csv | tee psauron_protein_${base}.log

    { echo "--> PSAURON spliced-CDS (full model) =="; cat psauron_cds_${base}.log
      echo "--> PSAURON protein (-p) =="            ; cat psauron_protein_${base}.log
    } > psauron_summary_${base}.txt
    """

    stub:
    def base = cds.simpleName
    """
    touch psauron_cds_${base}.csv psauron_protein_${base}.csv psauron_summary_${base}.txt
    """
}

process COMPLEASM {
    tag   { "${proteins.simpleName}_${lineage}" }
    label 'compleasm'
    conda 'bioconda::compleasm'
    publishDir "${params.outdir}/compleasm", mode: 'copy'

    input:
    path proteins
    val  lineage
    path library
    path runner

    output:
    path "compleasm_${proteins.simpleName}", type: 'dir'
    path "compleasm_summary_${proteins.simpleName}.txt", emit: summary

    script:
    def base = proteins.simpleName
    """
    python3 ${runner} protein -p ${proteins} -l ${lineage} -t ${task.cpus} \
        -o compleasm_${base} -L ${library}
    cp compleasm_${base}/summary.txt compleasm_summary_${base}.txt
    # Drop the per-marker hmmsearch files (~11k tiny files per run). They are
    # intermediates; summary.txt + full_table.tsv are the real outputs.
    # Keeping them blows the scratch inode quota across multiple species.
    rm -rf compleasm_${base}/*_hmmsearch_output
    """

    stub:
    def base = proteins.simpleName
    """
    mkdir -p compleasm_${base}
    touch compleasm_${base}/summary.txt compleasm_summary_${base}.txt
    """
}

process OMARK {
    tag   { proteins.simpleName }
    label 'omark'
    conda 'bioconda::omark bioconda::omamer'
    publishDir "${params.outdir}/omark", mode: 'copy'

    input:
    path proteins
    path db
    path ete_db
    path make_isoform

    output:
    path "omark_${proteins.simpleName}", type: 'dir'
    path "omark_summary_${proteins.simpleName}.txt", emit: summary

    script:
    def base = proteins.simpleName
    """
    omamer search --db ${db} --query ${proteins} --out ${base}.omamer
    mkdir -p omark_${base}
    # Group isoforms per gene so OMArk keeps only the best isoform and does not
    # count alternative transcripts as "Duplicated" completeness.
    python3 ${make_isoform} ${proteins} ${base}.splice
    # -e points OMArk at a prebuilt ete3 NCBI taxonomy DB so it never tries
    # to reach the internet (compute nodes have no outbound network).
    omark -f ${base}.omamer -d ${db} -e ${ete_db} -i ${base}.splice -o omark_${base}
    cat omark_${base}/*summary*.txt > omark_summary_${base}.txt
    """

    stub:
    def base = proteins.simpleName
    """
    mkdir -p omark_${base}
    touch omark_${base}/${base}_summary.txt omark_summary_${base}.txt
    """
}

process AGAT_STATS {
    tag   { gff.simpleName }
    label 'agat'
    conda 'bioconda::agat'
    publishDir "${params.outdir}/agat", mode: 'copy'

    input:
    path gff
    path genome

    output:
    path "agat_stats_${gff.simpleName}.txt",     emit: stats
    path "agat_premature_${gff.simpleName}.gff3", emit: premature

    script:
    def base = gff.simpleName
    """
    agat_sp_statistics.pl --gff ${gff} -o agat_stats_${base}.txt
    # internal / premature stop codons (needs genome to translate CDS)
    agat_sp_flag_premature_stop_codons.pl --gff ${gff} --fasta ${genome} \
        -o agat_premature_${base}.gff3
    """

    stub:
    def base = gff.simpleName
    """
    touch agat_stats_${base}.txt agat_premature_${base}.gff3
    """
}

process GFFUTILS {
    tag   { gff.simpleName }
    label 'python'
    conda 'conda-forge::python bioconda::gffutils'
    publishDir "${params.outdir}/gffutils", mode: 'copy'

    input:
    path gff
    path script

    output:
    path "gffutils_${gff.simpleName}.txt", emit: txt

    script:
    """
    python3 ${script} ${gff} > gffutils_${gff.simpleName}.txt
    """

    stub:
    """
    touch gffutils_${gff.simpleName}.txt
    """
}

process PLOT_GENELENGTH {
    tag   { gff.simpleName }
    label 'plot'
    conda 'conda-forge::python conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn'
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path gff
    path script

    output:
    path "genelength_${gff.simpleName}.png"

    script:
    """
    python3 ${script} ${gff} genelength_${gff.simpleName}.png
    """

    stub:
    """
    touch genelength_${gff.simpleName}.png
    """
}

process PLOT_GENEDENSITY {
    tag   { gff.simpleName }
    label 'plot'
    conda 'conda-forge::python conda-forge::pandas conda-forge::matplotlib'
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path gff
    path script

    output:
    path "genedensity_${gff.simpleName}.png"

    script:
    """
    python3 ${script} ${gff} genedensity_${gff.simpleName}.png
    """

    stub:
    """
    touch genedensity_${gff.simpleName}.png
    """
}

process PLOT_CPG {
    tag   { genome.simpleName }
    label 'plot_r'
    conda 'bioconda::bioconductor-karyoploter bioconda::bioconductor-regioner bioconda::bioconductor-biostrings bioconda::bioconductor-genomicranges'
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path genome
    path script

    output:
    path "cpg_${genome.simpleName}.png"
    path "cpg_${genome.simpleName}.tsv"

    script:
    """
    Rscript ${script} ${genome} cpg_${genome.simpleName}.png cpg_${genome.simpleName}.tsv
    """

    stub:
    """
    touch cpg_${genome.simpleName}.png cpg_${genome.simpleName}.tsv
    """
}

process COMPARE_REF {
    tag   { "${pred.simpleName}_vs_${ref.simpleName}" }
    label 'agat'
    conda 'bioconda::agat bioconda::gffcompare bioconda::aegean conda-forge::python conda-forge::pandas conda-forge::matplotlib conda-forge::seaborn'
    publishDir "${params.outdir}/compare_ref", mode: 'copy'

    input:
    path ref
    path pred
    path script

    output:
    path "compare_${pred.simpleName}", type: 'dir'
    path "compare_summary_${pred.simpleName}.txt", emit: summary

    script:
    def base = pred.simpleName
    """
    mkdir -p compare_${base}

    # 1. AGAT structural comparison
    agat_sp_compare_two_annotations.pl -gff1 ${ref} -gff2 ${pred} \
        -o compare_${base}/agat_compare || true

    # 2. gffcompare sensitivity / precision
    gffcompare -r ${ref} -o compare_${base}/gffcompare ${pred}

    # 3. parseval (normalise both first; write to NEW names, never over input)
    canon-gff3 ${ref}  > compare_${base}/ref.canon.gff3
    canon-gff3 ${pred} > compare_${base}/pred.canon.gff3
    parseval --summary --refrlabel reference --predlabel new \
        compare_${base}/ref.canon.gff3 compare_${base}/pred.canon.gff3 \
        > compare_${base}/parseval.txt || true

    # 4. overlaid gene-length distribution + stats
    python3 ${script} ${pred} ${ref} compare_${base}/compare

    cat compare_${base}/gffcompare*.stats compare_${base}/parseval.txt \
        compare_${base}/compare_stats.tsv > compare_summary_${base}.txt 2>/dev/null || \
        touch compare_summary_${base}.txt
    """

    stub:
    def base = pred.simpleName
    """
    mkdir -p compare_${base}
    touch compare_${base}/gffcompare.stats compare_summary_${base}.txt
    """
}

process SUMMARY_REPORT {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path reports

    output:
    path "SUMMARY_REPORT.txt"

    script:
    """
    {
      echo "# post-annotation evaluation summary"
      echo "# generated: \$(date)"
      for f in ${reports}; do
        echo
        echo "--------------------------------------------------------"
        echo "== \$f"
        echo "--------------------------------------------------------"
        cat \$f
      done
    } > SUMMARY_REPORT.txt
    """

    stub:
    """
    touch SUMMARY_REPORT.txt
    """
}
