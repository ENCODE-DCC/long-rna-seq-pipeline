<!-- dx-header -->
# ENCODE small-rna-seq-pipeline (DNAnexus Platform App)

This folder contains the dnanexis applets used in the ENCODE small-rna-seq pipeline. All applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

## Applets (aka steps):
- small-rna-prep-star - Takes a (gender specific) genome reference gzipped fasta file and produces a STAR genome 
                        index tar.gz file. This 'index' file needs to be produced one and will be used as input for 
                        every subsequent small-rna-seq run for the same genome reference (e.g. GRCh37/hg19 female). 
- small-rna-align     - Takes a (single-end) gzipped fastq file and the STAR genome index tar.gz file (created by 
                        the '`small-rna-prep-star`' applet).  This step produces a bam file of the aligned reads,
                        and a gene quantification file. 
- small-rna-signals   - Takes the bam file output from '`small-rna-align`' and a chromosome name/length file to produce 
                        four bigWig "signal' files for easy display in a genome browser.  The alignment signals are 
                        filtered by +/- strand and uniquely mapped vs. all mapped reads.
- small-rna-mad-qc    - Takes two gene quantification files produced from '`small-rna-align`' and calculates the Mean 
                        Absolute Deviation and correlations. This step produces a plot (png) file and some QC metric values.
                     
---------
## Flow
```
INPUTS:  genome.fasta.gz           reads.fq.gz            small-rna.bam(b)  2*quantification.cvs(c)
                                   star-index.tgz(a)      chrom.sizes                 |
                |                  chrom.sizes                  |                     |
                |                       |                       |                     |
                V                       V                       V                     V
STEPS:   small-rna-prep-star ====> small-rna-align ====> small-rna-signals ===> small-rna-mad-qc
                |                       |                       |                     |
                V                       V                       V                     V
OUTPUTS: star-index.tgz(a)         small-rna.bam(b)      unique-plus.bw            plot.png
                                   quantification.cvs(c) unique-minus.bw
                                                         all-plus.bw
                                                         all-minus.bw
```
---------
## Pipeline Restrictions for STAR based ENCODE small RNA-seq processing pipeline, <200bp

* The read length should be a minimum of 50 base pairs.
* Sequencing should be single-ended.
* All Illumina platforms are supported for use in the uniform pipeline; colorspace (SOLiD) are not supported.
* Barcodes and spike-in sequences, if present in the fastq, must be indicated.
* Library insert size range must be indicated. 
* Alignment files are mapped to either the GRCh38 or mm10 (https://www.encodeproject.org/references/ENCSR425FOI/) sequences.
* Gene and transcript quantification files are annotated to either GENCODE V24 or M4 (https://www.encodeproject.org/references/ENCSR884DHJ/)
                                                         
