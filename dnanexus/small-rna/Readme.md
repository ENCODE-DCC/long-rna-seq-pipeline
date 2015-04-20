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
                        the '`small-rna-prep-star`' applet).  This step produces a bam file of the aligned reads. 
- small-rna-signals   - Takes the bam file output from '`small-rna-align`' and a chromosome name/length file to produce 
                        four bigWig "signal' files for easy display in a genome browser.  The alignment signals are 
                        filtered by +/- strand and uniquely mapped vs. all mapped reads.
                     
---------
## Flow
```
INPUTS:  genome.fasta.gz           reads.fq.gz            rampage.bam(a)
                                   star-index.tgz(a)      chrom.sizes
                |                  chrom.sizes                  |
                |                       |                       |
                V                       V                       V
STEPS:   small-rna-prep-star ====> small-rna-align ====> small-rna-signals
                |                       |                       |
                V                       V                       V
OUTPUTS: star-index.tgz(a)         small-rna.bam(b)      unique-plus.bw
                                                         unique-minus.bw
                                                         all-plus.bw
                                                         all-minus.bw
```
                                                         
