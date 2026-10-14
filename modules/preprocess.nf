/*
 * Per sample: primer trimming + orientation, removal of reads with a second
 * full primer (concatemers), extraction of the full ITS by cutting at the
 * conserved ITS4 site, and read renaming as <read_id>;sample=<sample>.
 * With its_end_site = 'none' (e.g. COI or other amplicons with no conserved
 * internal site) reads are not cut, only filtered by min/max_its_len.
 */
process PREPROCESS {
    tag "${sample}"
    label 'process_medium'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    tuple val(sample), path(reads, stageAs: 'input/*')

    output:
    tuple val(sample), path("${sample}_ITS.fastq.gz"), emit: its
    path "${sample}_read_counts.tsv",                  emit: counts
    path "${sample}_ITS_stats.tsv",                    emit: stats
    path "${sample}_cutadapt_*.log",                   emit: logs
    path "versions.yml",                               emit: versions

    script:
    def rc_fwd   = revcomp(params.primer_fwd)
    def rc_rev   = revcomp(params.primer_rev)
    def rev_mode = (params.require_rev_primer as String) == 'true' ? 'required' : 'optional'
    def ov_fwd   = params.primer_fwd.length() - 4
    def ov_rev   = params.primer_rev.length() - 4
    def cut_its  = (params.its_end_site as String).toLowerCase() != 'none'
    def ov_its   = cut_its ? params.its_end_site.length() - 5 : 0
    """
    count() { seqkit stats -T "\$1" | awk -F'\\t' 'NR==2{print \$4}'; }
    printf "sample\\tstep\\treads\\n" > ${sample}_read_counts.tsv

    # all input files (one file, or the split files of a barcode folder) -> one file
    seqkit seq input/* -o raw.fastq.gz
    printf "%s\\traw\\t%s\\n" ${sample} \$(count raw.fastq.gz) >> ${sample}_read_counts.tsv

    # 1. primers: linked F...rc(R); --rc orients every read F->R; anything outside
    #    the primers (adapters, barcodes) is removed
    cutadapt -j ${task.cpus} -e ${params.primer_error} --rc \\
        -g "${params.primer_fwd};min_overlap=${ov_fwd}...${rc_rev};min_overlap=${ov_rev};${rev_mode}" \\
        --discard-untrimmed -m ${params.min_amplicon_len} -M ${params.max_amplicon_len} \\
        -o primer_trimmed.fastq.gz raw.fastq.gz > ${sample}_cutadapt_primers.log
    printf "%s\\tprimer_trimmed\\t%s\\n" ${sample} \$(count primer_trimmed.fastq.gz) >> ${sample}_read_counts.tsv

    # 2. a second full primer inside an amplicon indicates a concatemer/chimera
    cutadapt -j ${task.cpus} -e ${params.primer_error} \\
        -b "F=${params.primer_fwd};min_overlap=${params.primer_fwd.length()}" \\
        -b "rcF=${rc_fwd};min_overlap=${params.primer_fwd.length()}" \\
        -b "R=${params.primer_rev};min_overlap=${params.primer_rev.length()}" \\
        -b "rcR=${rc_rev};min_overlap=${params.primer_rev.length()}" \\
        --discard-trimmed -o no_internal_primers.fastq.gz primer_trimmed.fastq.gz \\
        > ${sample}_cutadapt_internal.log
    printf "%s\\tno_internal_primers\\t%s\\n" ${sample} \$(count no_internal_primers.fastq.gz) >> ${sample}_read_counts.tsv

    # 3. full ITS: cut at the ITS4 site (~40 bp into the LSU); with its_end_site
    #    'none', reads are kept whole and only filtered by length
    if [ "${cut_its}" = true ]; then
        cutadapt -j ${task.cpus} -e 0.2 \\
            -a "${params.its_end_site};min_overlap=${ov_its}" \\
            --discard-untrimmed -m ${params.min_its_len} -M ${params.max_its_len} \\
            -o its_tmp.fastq.gz no_internal_primers.fastq.gz > ${sample}_cutadapt_its.log
    else
        seqkit seq -m ${params.min_its_len} -M ${params.max_its_len} \\
            -o its_tmp.fastq.gz no_internal_primers.fastq.gz
        echo "its_end_site = none: reads not cut; kept ${params.min_its_len}-${params.max_its_len} bp" > ${sample}_cutadapt_its.log
    fi

    # 4. rename reads with the sample
    seqkit replace -p '\\s.*\$' -r ";sample=${sample}" -o ${sample}_ITS.fastq.gz its_tmp.fastq.gz
    printf "%s\\tfull_ITS\\t%s\\n" ${sample} \$(count ${sample}_ITS.fastq.gz) >> ${sample}_read_counts.tsv
    seqkit stats -a -T ${sample}_ITS.fastq.gz > ${sample}_ITS_stats.tsv
    rm raw.fastq.gz primer_trimmed.fastq.gz no_internal_primers.fastq.gz its_tmp.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}

/*
 * Optional: SINTAX classification of every full-ITS read.
 */
process PER_READ_SINTAX {
    tag "${sample}"
    label 'process_high'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    tuple val(sample), path(reads)
    path db

    output:
    path "${sample}_per_read_sintax.tsv.gz", emit: tsv
    path "versions.yml",                     emit: versions

    script:
    """
    vsearch --sintax ${reads} --db ${db} --sintax_cutoff ${params.sintax_cutoff} \\
        --strand plus --randseed ${params.sintax_seed} --threads ${task.cpus} --tabbedout ${sample}_per_read_sintax.tsv
    gzip ${sample}_per_read_sintax.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
    END_VERSIONS
    """
}

def revcomp(String seq) {
    return seq.reverse().tr('ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn', 'TGCAAWSKMYRVHDBNtgcaawskmyrvhdbn')
}
