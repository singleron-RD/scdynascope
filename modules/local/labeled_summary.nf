process LABELED_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::scanpy==1.10.2 conda-forge::numpy==1.26.4'
    container "docker://gcfntnu/scanpy:1.10.1-dev"

    input:
    tuple val(meta), path(labeled_dir), path(filtered_dir)
    val min_cells
    val min_genes


    output:
    tuple val(meta), path("*.json"),        emit: json
    //path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    labeled_summary.py \\
        --outdir ./ \\
        --sample ${meta.id} \\
        --filtered_matrix ${filtered_dir} \\
        --labeled_matrix ${labeled_dir} \\
        --min_cells ${min_cells} \\
        --min_genes ${min_genes}
    """

}