/*
 * Consensus polishing of every OTU centroid with its member reads (medaka).
 * Equivalent to running medaka_consensus once per OTU, but batched: each OTU's
 * reads are aligned only to its own centroid (minimap2 with the model's
 * alignment settings), the alignments are merged into one BAM, and medaka
 * inference + sequence run once, so the model is loaded only once.
 * OTUs medaka cannot polish keep their centroid (polished=no).
 */
process POLISH {
    label 'process_gpu'
    conda     params.conda_medaka ?: "${moduleDir}/../envs/medaka.yml"
    container params.container_medaka

    input:
    path clusters
    path otu_map

    output:
    path "consensus.fasta",        emit: consensus
    path "medaka_*.log"
    path "versions.yml",           emit: versions

    script:
    """
    python3 -c "import torch; print('torch', torch.__version__, '| CUDA GPU available:', torch.cuda.is_available())" \\
        > medaka_device.log 2>&1 || true

    OTUS=\$(awk -F'\\t' 'NR>1{print \$1}' ${otu_map})
    for otu in \$OTUS; do cat ${clusters}/\${otu}_centroid.fasta; done > drafts.fasta
    samtools faidx drafts.fasta

    # medaka's model-specific scoring in mini_align terms (-M -S -O -E) -> minimap2 (-A -B -O -E)
    SCORING=\$(medaka tools get_alignment_params --model ${params.medaka_model} 2>/dev/null | tail -1 \\
        | sed -e 's/-M /-A /' -e 's/-S /-B /')
    {
        awk -F'\\t' '{print "@SQ\\tSN:"\$1"\\tLN:"\$2}' drafts.fasta.fai
        for otu in \$OTUS; do
            minimap2 -x map-ont --secondary=no -L \$SCORING -t ${task.cpus} -a \\
                ${clusters}/\${otu}_centroid.fasta ${clusters}/\${otu}.fastq 2>/dev/null \\
                | samtools view -F 2308 -
        done
    } | samtools sort -@ ${task.cpus} -o calls_to_draft.bam -
    samtools index calls_to_draft.bam

    medaka inference calls_to_draft.bam consensus_probs.hdf \\
        --model ${params.medaka_model} --batch_size 100 --threads ${task.cpus} > medaka_inference.log 2>&1
    medaka sequence consensus_probs.hdf drafts.fasta medaka_consensus.fasta \\
        --threads ${task.cpus} > medaka_sequence.log 2>&1

    collect_consensus.py medaka_consensus.fasta drafts.fasta ${otu_map} consensus.fasta | tee -a medaka_device.log
    rm -f consensus_probs.hdf calls_to_draft.bam*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version | sed 's/medaka //')
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
