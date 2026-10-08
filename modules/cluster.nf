/*
 * All samples pooled: sort reads by mean quality (drop < min_read_q), cluster in
 * decreasing quality order so the best reads become centroids, and write the
 * member reads of every cluster with >= min_cluster_size reads for polishing.
 */
process CLUSTER {
    label 'process_high'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path its_reads

    output:
    path "clusters",                emit: clusters
    path "otu_map.tsv",             emit: otu_map
    path "otu_counts.tsv",          emit: counts
    path "cluster_summary.tsv",     emit: summary
    path "clusters.uc.gz"
    path "centroids.fasta.gz"
    path "read_quality.tsv.gz"
    path "vsearch_cluster.log"
    path "versions.yml",            emit: versions

    script:
    """
    seqkit seq ${its_reads} > pooled.fastq

    sort_reads_by_quality.py -i pooled.fastq -q ${params.min_read_q} \\
        -o pooled_qualsorted.fastq -t read_quality.tsv
    rm pooled.fastq

    # reads are already oriented, so --strand plus
    vsearch --cluster_smallmem pooled_qualsorted.fastq --usersort \\
        --id ${params.cluster_id} --strand plus --threads ${task.cpus} \\
        --centroids centroids.fasta --uc clusters.uc 2> vsearch_cluster.log

    make_cluster_files.py -c clusters.uc -f pooled_qualsorted.fastq -o . \\
        -s ${params.min_cluster_size} -m ${params.max_reads_per_cluster}

    rm pooled_qualsorted.fastq*
    gzip clusters.uc centroids.fasta read_quality.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
