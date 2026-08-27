## Evaluation of eviann genome annotation
---

## NextFlow Evaluation Pipeline
#### Evaluation of eviann genome annotation
---

**_Tools_**
- `PSAUSEON` -> `compleasm` -> `OMarK` ->
- `AGAT` -> `GFFUTILS`
- gene density -> gene length -> cpg islands
  => **SUMMARY REPORT**

---

##### from `main.nf` -> **_`params & workflow`_**

```groovy
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
```
