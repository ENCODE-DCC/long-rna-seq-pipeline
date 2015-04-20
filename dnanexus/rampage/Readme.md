<!-- dx-header -->
# ENCODE rampage-rna-seq-pipeline (DNAnexus Platform App)

This folder contains the dnanexis applets used in the ENCODE rampage-rna-seq pipeline. All applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

## Applets (aka steps):
- rampage-align-pe - Takes a pair of (paired-end) gzipped fastq files and a STAR genome index tar.gz file (identical 
                     to and created by the long-rna-seq-pipeline's '`prep-star`' applet).  This step produces a bam file 
                     of the aligned reads. 
- rampage-signals  - Takes the bam file output from '`ramapge-align-pe`' and a chromosome name/length file to produce 
                     four bigWig 'signal' files for easy display in a genome browser.  The 5' TSS alignment signals
                     are filtered by +/- strand and uniquely mapped vs. all mapped reads.
- rampage-peaks    - Takes the bam file output from '`ramapge-align-pe`' and a second 'control' genome-aligned bam 
                     (produced by the long-rna-seq-pipeline's '`align-star-pe`' applet).  An annotation gzipped gtf 
                     file (such as gencode.v19.annotation.gtf.gz) and a chromosome name/length file is also required. 
                     This peak calling step produces three files all with essentially the same data in different 
                     formats: gff, bed, and bigBed for easy display in a genome browser.  
- rampage idr      - Takes the bed output from two separate runs of '`rampage-peaks`' and performs an "Irreproducible 
                     Discover Rate" analysis comparing the 2 sets of peaks.  The input peaks are typically from a 
                     pair of replicates of the same experiment allowing validation of the experiment methods, as well
                     as a means of reducing noise in the final results.  This step also needs a chromosome name/size 
                     file and produces a bed and a bigWig file, each containing the superset of peaks and the idr 
                     statistics on the overlapping set of peaks.  This applet also produces a plot in png format.
                     
---------
## Flow
```
INPUTS:  read1.fq.gz            rampage.bam(a)    rampage.bam(a)     peaks.bed(b) from replicate 1
         read2.fq.gz            chrom.sizes       long-rna.bed       peaks.bed(b) from replicate 2
         star-index.tar.gz          |             chrom.sizes        chrom.sizes
            |                       |             annotation.gtf.gz      |
            |                       |                 |                  |
            V                       V                 V                  V
STEPS:   rampage-align-pe ====> rampage-signals,  rampage-peaks ===> rampage-idr
            |                       |                 |                  |
            V                       V                 V                  V
OUTPUTS: rampage.bam(a)         unique-plus.bw    peaks.gtf          peaks-idr.bed
                                unique-minus.bw   peaks.bed(b)       peaks-idr.bb
                                all-plus.bw       peaks.bb           idr-plot.png
                                all-minus.bw
```

