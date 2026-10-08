/*
 * Merge near-identical polished consensus OTUs, annotate and classify them, and
 * build the final tables.
 */

process MERGE_OTUS {
    label 'process_low'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus
    path counts

    output:
    path "merged_consensus.fasta",   emit: consensus
    path "merged_otu_counts.tsv",    emit: counts
    path "otu_membership.tsv",       emit: membership, optional: true
    path "versions.yml",             emit: versions

    script:
    if (params.merge_id as double >= 1.0)
        """
        cp ${consensus} merged_consensus.fasta
        cp ${counts} merged_otu_counts.tsv
        echo '"${task.process}": {}' > versions.yml
        """
    else
        """
        # Read-level clustering splits one template into many OTUs (reads with more
        # error than 1 - cluster_id vs. the centroid start new clusters); after
        # polishing these are near-identical and are merged by abundance here.
        vsearch --cluster_size ${consensus} --sizein --id ${params.merge_id} --strand plus \\
            --threads ${task.cpus} --uc merge_consensus.uc 2> merge_consensus.log
        merge_otus.py --consensus ${consensus} --counts ${counts} --uc merge_consensus.uc --outdir merged
        mv merged/consensus.fasta merged_consensus.fasta
        mv merged/otu_counts.tsv merged_otu_counts.tsv
        mv merged/otu_membership.tsv otu_membership.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
        END_VERSIONS
        """
}

process ITSX {
    label 'process_medium'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus

    output:
    path "itsx.positions.txt",  emit: positions
    path "itsx.*"
    path "versions.yml",        emit: versions

    script:
    """
    ITSx -i ${consensus} -o itsx -t ${params.itsx_taxa} --cpu ${task.cpus} \\
        --save_regions all --graphical F --preserve T --complement F > itsx.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ITSx: \$(ITSx --help 2>&1 | grep -m1 -i '^Version' | sed 's/Version: //')
        hmmer: \$(hmmscan -h | grep -m1 HMMER | sed 's/# HMMER //; s/ .*//')
    END_VERSIONS
    """
}

process CHIMERA {
    label 'process_low'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus

    output:
    path "uchime3.tsv",     emit: tsv
    path "versions.yml",    emit: versions

    script:
    """
    vsearch --uchime3_denovo ${consensus} --uchimeout uchime3.tsv 2> uchime3.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
    END_VERSIONS
    """
}

process TOPHIT {
    label 'process_high'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus
    path db

    output:
    path "top_hit.tsv",     emit: tsv
    path "versions.yml",    emit: versions

    script:
    """
    # best global-alignment hit with >= tax_query_cov query coverage; several
    # candidates are examined so the best (not the first acceptable) hit is kept
    vsearch --usearch_global ${consensus} --db ${db} \\
        --id 0.5 --query_cov ${params.tax_query_cov} --strand plus \\
        --maxaccepts 50 --maxrejects 256 --maxhits 1 --top_hits_only \\
        --threads ${task.cpus} --output_no_hits \\
        --userout top_hit.tsv --userfields query+target+id+alnlen+qcov+tcov+mism+opens 2> top_hit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
    END_VERSIONS
    """
}

process SINTAX {
    label 'process_high'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus
    path db

    output:
    path "sintax.tsv",      emit: tsv
    path "versions.yml",    emit: versions

    script:
    """
    vsearch --sintax ${consensus} --db ${db} --sintax_cutoff ${params.sintax_cutoff} \\
        --strand plus --threads ${task.cpus} --tabbedout sintax.tsv 2> sintax.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/^vsearch v//; s/,.*//')
    END_VERSIONS
    """
}

process OTU_TABLE {
    label 'process_low'
    conda     params.conda_main ?: "${moduleDir}/../envs/main.yml"
    container params.container_main

    input:
    path consensus
    path counts
    path top_hit
    path sintax
    path itsx_positions
    path uchime
    path membership
    path read_counts
    path cluster_summary

    output:
    path "OTU_table.tsv"
    path "OTU_info.tsv"
    path "OTU_consensus.fasta"
    path "OTU_merge_membership.tsv", optional: true
    path "read_tracking.tsv"
    path "versions.yml",             emit: versions

    script:
    """
    make_otu_table.py --consensus ${consensus} --counts ${counts} \\
        --top_hit ${top_hit} --sintax ${sintax} --itsx_positions ${itsx_positions} \\
        --uchime ${uchime} --outdir .
    if [ -s "${membership}" ] && [ "${membership}" != "NO_FILE" ]; then
        cp ${membership} OTU_merge_membership.tsv
    fi
    read_tracking.py --cluster_summary ${cluster_summary} --counts ${read_counts} > read_tracking.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
