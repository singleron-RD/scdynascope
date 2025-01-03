


- [Main Output](#main-output)
- [modules](#modules)
  - [filter\_gtf](#filter_gtf)
  - [star\_genome](#star_genome)
  - [protocol\_cmd](#protocol_cmd)
  - [starsolo](#starsolo)
  - [cell\_calling](#cell_calling)
  - [starsolo\_summary](#starsolo_summary)
  - [split\_bam](#split_bam)
  - [conversion](#conversion)
  - [conversion\_merge](#conversion_merge)
  - [substitution](#substitution)
  - [labeled](#labeled)
  - [labeled\_summary](#labeled_summary)
  - [multiqc-sgr](#multiqc-sgr)
  - [pipeline\_info](#pipeline_info)
  - [subsample(Optional)](#subsampleoptional)
  - [fastqc(Optional)](#fastqcoptional)

# Main Output

- `multiqc/multiqc_report.html` HTML report containing QC metrics.
- `cell_calling/{sample}.matrix/filtered` Gene expression matrix file contains only cell barcodes. This file should be used as input to downstream analysis tools such as Seurat and Scanpy.
- `labeled/{sample}.matrix/labeled` and `labeled/{sample}.matrix/unlabeled` Labeled and unlabeled expression matrix. Gene expression matrix file contains only cell barcodes.
- `starsolo/{sample}.Aligned.sortedByCoord.out.bam` Bam file contains coordinate-sorted reads aligned to the genome.
- `substitution/{sample}.substitution.tsv` The overall substitution rates for each conversion type in reads.

# modules

## filter_gtf

This module has the same functionality as [`cellranger mkgtf`](https://kb.10xgenomics.com/hc/en-us/articles/360002541171-What-criteria-should-I-use-with-the-mkgtf-tool-when-making-a-custom-reference-for-Cell-Ranger)

> GTF files can contain entries for non-polyA transcripts that overlap with protein-coding gene models. These entries can cause reads to be flagged as mapped to multiple genes (multi-mapped) because of the overlapping annotations. In the case where reads are flagged as multi-mapped, they are not counted.

> We recommend filtering the GTF file so that it contains only gene categories of interest by using the cellranger mkgtf tool. Which genes to filter depends on your research question.

The filtering criteria is controlled by the argument `--keep_attributes`. The default value of this argument is the same as the [reference used by cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_3.0.0)

> [!NOTE]
> gtf files from [genecode](https://www.gencodegenes.org/) use `gene_type` instead of `gene_biotype`.
>
> ```
> --keep_attributes "gene_type=protein_coding,lncRNA..."
> ```

**Output files**

- `*.filtered.gtf` GTF file after filtering.
- `gtf_filter.log` log file containing number of lines filtered in the original gtf file.

## star_genome

Generate STAR genome index. Detailed documents can be found in the [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

> [!TIP]
> Once you have the indices from a workflow run you should save them somewhere central and reuse them in subsequent runs using custom config files or command line parameters.

**Output files**

- `{genome_name}/` STAR genome index folder.

## protocol_cmd

Automatically detect [GEXSCOPE protocol](../assets/protocols.json) from R1 reads and generate STARSolo command-line arguments accordingly.

**Output files**

- `{sample}.protocol.txt` Detected protocol.
- `{sample}.starsolo_cmd.txt` STARSolo command-line arguments.

## starsolo

Descriptions of parameters and files can be found in [STARSolo documents](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).
When you have questions, [STAR’s github issue](https://github.com/alexdobin/STAR/issues) is also a great place to find answers and help.

> [!NOTE]
> The command line arguments in this STARsolo documentation may not be up to date. For the latest STARSolo arguments, please refer to The STAR Manual.

**Output files**

- `{sample}.Aligned.sortedByCoord.out.bam` Bam file contains coordinate-sorted reads aligned to the genome.

## cell_calling
Run the cell filtering algorithm on the previously generated raw matrix.

**Output files**

- `{sample}.matrix` Raw and filtering matrix.

## starsolo_summary

Extract data for visualization from starsolo result files.

## split_bam

Split bam file for downstream analysis in parallel.

## conversion

Get conversion for each read and add tags. 

**Output files**

- `{sample}.PosTag.bam` BAM file with conversion tags added.
- `{sample}.PosTag.csv` A list of sites of the specified conversion type.
- `{sample}.snp.csv` Candidate snp sites obtained according to the screening conditions.

## substitution

Computes the overall conversion rates in reads and plots a barplot.

**Output files**

- `{sample}.TC_substitution.tsv` Substitution rate of labeled and background sites.

## labeled

Quantify unlabeled and labeled RNA. (Labeled samples only)

**Output files**

- `{sample}.matrix/labeled` labeled matrix
- `{sample}.matrix/unlabeled` unlabeled matrix

## labeled_summary

Boxplots for TOR rates distribution of cells or genes.

## multiqc-sgr

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

[multiqc-sgr](https://pypi.org/project/multiqc-sgr/) adds some modules on this basis to facilitate the visualization of single cell-related data.

**Output files**

- `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
- `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
- `multiqc_plots/`: directory containing static images from the report in various formats.

## pipeline_info

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files**

- Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
- Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
- Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
- Parameters used by the pipeline run: `params.json`.

## subsample(Optional)
Randomly extract 10%, 20%, ..., 80%, 90% of the reads in the bam file and calculate the corresponding saturation and median genes respectively.

**Output files**

Saturation and median genes plots are added to the multiqc report.


## fastqc(Optional)

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output files**

- `*_fastqc.html`: FastQC report containing quality metrics.
- `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
