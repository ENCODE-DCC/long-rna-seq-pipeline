# long-rna-seq-pipeline
STAR based ENCODE Long RNA-Seq processing pipeline: **Pipeline Restrictions**

* The read length should be a minimum of 50 base pairs.
* Sequencing may be paired- or single-end, as long as sequencing type is specified and read pairs are indicated.
* All Illumina platforms are supported for use in the uniform pipeline; colorspace (SOLiD) are not supported.
* Barcodes, if present in fastq, must be indicated in the metadata.
* ERCC spike-ins should be used in library preparation with the concentrations indicated in the metadata. 
* Library insert size range must be indicated. 
* Alignment files are mapped to either the GRCh38 or mm10 sequences.
* Gene and transcript quantification files are annotated to either GENCODE V24 or M4.
