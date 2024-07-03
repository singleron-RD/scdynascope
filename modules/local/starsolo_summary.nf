process STARSOLO_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    input:
<<<<<<< HEAD
    tuple val(meta), path(read_stats)
    path summary
    path barcodes
=======
    tuple val(meta), path(read_stats), path(summary), path(filtered_matrix)
>>>>>>> upstream/master

    output:
    tuple val(meta), path("*.json"), emit: json

    script:

    """
    starsolo_summary.py \\
        --read_stats ${read_stats} \\
<<<<<<< HEAD
        --barcodes ${barcodes} \\
=======
        --filtered_matrix ${filtered_matrix} \\
>>>>>>> upstream/master
        --summary ${summary} \\
        --sample ${meta.id}
    """
}