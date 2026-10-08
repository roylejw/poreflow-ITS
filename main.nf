#!/usr/bin/env nextflow
/*
 * poreflow-ITS - full-length ITS metabarcoding of Oxford Nanopore amplicon reads
 *
 *   primer trimming + orientation -> full-ITS extraction (ITS4 site)
 *   -> quality-ordered clustering -> medaka consensus per OTU
 *   -> merging of near-identical consensus OTUs -> ITSx, chimera check,
 *      top hit + SINTAX taxonomy -> OTU table
 *
 * Adapted from Aline-Git/18S_metabarcoding_pipeline (GPL-3.0).
 */

include { PREPROCESS; PER_READ_SINTAX }                            from './modules/preprocess'
include { CLUSTER }                                                from './modules/cluster'
include { POLISH }                                                 from './modules/polish'
include { MERGE_OTUS; ITSX; CHIMERA; TOPHIT; SINTAX; OTU_TABLE }   from './modules/taxonomy'

workflow {
    main:
    if ((params.help as String) == 'true') {
        log.info(helpMessage())
    }
    else {
        validateParams()
        log.info(summary())

        def db = file(params.db, checkIfExists: true)
        samples = channel.fromList(parseInput(params.input as String))

        PREPROCESS(samples)
        if ((params.per_read_tax as String) == 'true') {
            PER_READ_SINTAX(PREPROCESS.out.its, db)
        }

        CLUSTER(PREPROCESS.out.its.map { _sample, reads -> reads }.collect())
        POLISH(CLUSTER.out.clusters, CLUSTER.out.otu_map)
        MERGE_OTUS(POLISH.out.consensus, CLUSTER.out.counts)

        ITSX(MERGE_OTUS.out.consensus)
        CHIMERA(MERGE_OTUS.out.consensus)
        TOPHIT(MERGE_OTUS.out.consensus, db)
        SINTAX(MERGE_OTUS.out.consensus, db)

        OTU_TABLE(
            MERGE_OTUS.out.consensus,
            MERGE_OTUS.out.counts,
            TOPHIT.out.tsv,
            SINTAX.out.tsv,
            ITSX.out.positions,
            CHIMERA.out.tsv,
            MERGE_OTUS.out.membership.ifEmpty(file("${projectDir}/assets/NO_FILE")),
            PREPROCESS.out.counts.collect(),
            CLUSTER.out.summary
        )

        PREPROCESS.out.versions.first()
            .mix(CLUSTER.out.versions, POLISH.out.versions, MERGE_OTUS.out.versions, ITSX.out.versions,
                 CHIMERA.out.versions, TOPHIT.out.versions, SINTAX.out.versions, OTU_TABLE.out.versions)
            .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
    }
}

/*
 * Input: a samplesheet (CSV or TSV with columns sample,reads; reads = a fastq(.gz)
 * file or a folder of them, e.g. fastq_pass/barcode01) or a folder containing one
 * fastq(.gz) per sample and/or one sub-folder per sample.
 */
def parseInput(String input) {
    def fastq = ~/.+\.(fastq|fq)(\.gz)?$/
    def target = file(input, checkIfExists: true)
    def rows = []
    if (target.isDirectory()) {
        target.listFiles().sort { x -> x.name }.each { f ->
            if (f.isDirectory()) {
                def reads = f.listFiles().findAll { x -> x.name ==~ fastq }.sort { x -> x.name }
                if (reads) rows << [f.name, reads]
            }
            else if (f.name ==~ fastq) {
                rows << [f.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, ''), [f]]
            }
        }
    }
    else {
        target.readLines().eachWithIndex { line, i ->
            def fields = line.trim().split(/[,\t]/)*.trim()
            if (!line.trim() || line.startsWith('#') || (i == 0 && fields[0] == 'sample')) return
            if (fields.size() < 2) error("Samplesheet line ${i + 1}: expected 'sample,reads' but got: ${line}")
            def path = (fields[1].startsWith('/') || fields[1] ==~ /^[a-zA-Z0-9]+:\/\/.+/) ? file(fields[1]) : target.parent.resolve(fields[1])
            if (!path.exists()) error("Samplesheet line ${i + 1}: reads not found: ${path}")
            def reads = path.isDirectory() ? path.listFiles().findAll { x -> x.name ==~ fastq }.sort { x -> x.name } : [path]
            if (!reads) error("Samplesheet line ${i + 1}: no fastq files in ${path}")
            rows << [fields[0], reads]
        }
    }
    if (!rows) error("No samples found in ${input}")
    def names = rows.collect { r -> r[0] }
    def dups = names.findAll { n -> names.count(n) > 1 }.unique()
    if (dups) error("Duplicate sample names: ${dups.join(', ')}")
    def bad = names.findAll { n -> !(n ==~ /[A-Za-z0-9._-]+/) }
    if (bad) error("Sample names may only contain letters, digits, '.', '_' and '-': ${bad.join(', ')}")
    return rows
}

def validateParams() {
    if (!params.input) error("Please provide --input (samplesheet or folder of reads). See --help.")
    if (!params.db)    error("Please provide --db (SINTAX-formatted reference FASTA). See --help.")
    [cluster_id: params.cluster_id, merge_id: params.merge_id, tax_query_cov: params.tax_query_cov,
     sintax_cutoff: params.sintax_cutoff, primer_error: params.primer_error].each { k, v ->
        if (!(v as double >= 0 && v as double <= 1)) error("--${k} must be between 0 and 1 (got ${v})")
    }
    [primer_fwd: params.primer_fwd, primer_rev: params.primer_rev, its_end_site: params.its_end_site].each { k, v ->
        if (!(v ==~ /[ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]+/)) error("--${k} is not a DNA/IUPAC sequence: ${v}")
    }
}

def summary() {
    return """
    poreflow-ITS v${workflow.manifest.version}
    ==========================================================
    input            : ${params.input}
    reference db     : ${params.db}
    outdir           : ${params.outdir}
    primers          : ${params.primer_fwd} / ${params.primer_rev}
    ITS end site     : ${params.its_end_site}
    cluster / merge  : ${params.cluster_id} / ${params.merge_id} (min cluster size ${params.min_cluster_size})
    medaka model     : ${params.medaka_model}${(params.gpu as String) == 'true' ? ' (GPU)' : ''}
    profile          : ${workflow.profile}
    ==========================================================
    """.stripIndent()
}

def helpMessage() {
    return """
    poreflow-ITS v${workflow.manifest.version} - full-length ITS metabarcoding of Oxford Nanopore amplicon reads

    Usage:
      nextflow run main.nf -profile <executor>,<software> --input <samplesheet|folder> --db <reference.fasta> [options]
      nextflow run main.nf -profile <executor>,<software> -params-file my_params.yaml

    Profiles (combine one of each):
      executor : local, slurm, pbs, sge, lsf      (HPC profiles run every task in node-local \$TMPDIR)
      software : conda, mamba, docker, singularity, apptainer, wave
      test     : small bundled test run (see README)

    Required:
      --input      samplesheet (CSV/TSV with columns sample,reads) or a folder of reads
      --db         SINTAX-formatted reference FASTA (e.g. EUKARYOME or UNITE SINTAX release)

    Common options (defaults in brackets; full list in nextflow.config / params.example.yaml):
      --outdir               [${params.outdir}]
      --primer_fwd           [${params.primer_fwd}]
      --primer_rev           [${params.primer_rev}]
      --its_end_site         [${params.its_end_site}]
      --cluster_id           [${params.cluster_id}]
      --merge_id             [${params.merge_id}]   (1 = no merging of consensus OTUs)
      --min_cluster_size     [${params.min_cluster_size}]
      --medaka_model         [${params.medaka_model}]
      --gpu                  [${params.gpu}]   with --gpu_cluster_options e.g. '--gpus-per-node=1'
      --max_cpus / --max_memory / --max_time   [${params.max_cpus} / ${params.max_memory} / ${params.max_time}]
    """.stripIndent()
}
