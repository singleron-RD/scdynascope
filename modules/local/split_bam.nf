process SPLIT_BAM {
    tag "$meta.id"
    label 'process_single_high'

    conda "bioconda::picard=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam_sorted), path(log_final)
    val split_n_reads

    output:
    tuple val(meta), path("${meta.id}_split_bams/*.bam"), emit: split_bam_chunks
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def bam_chunks_dir = "${meta.id}_split_bams"
    def prefix = "${meta.id}"

    """
    total_reads="\$(grep 'Number of input reads |' ${log_final} | awk '{print \$6}')"
    total_reads=\$(( \${total_reads} > ${split_n_reads} ? \${total_reads} : ${split_n_reads} ))

    mkdir ${bam_chunks_dir}

    picard \\
        SplitSamByNumberOfReads \\
        I=${bam_sorted} \\
        O=${bam_chunks_dir} \\
        OUT_PREFIX=$prefix \\
        SPLIT_TO_N_READS=${split_n_reads} \\
        TOTAL_READS_IN_INPUT=\$total_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SplitSamByNumberOfReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
