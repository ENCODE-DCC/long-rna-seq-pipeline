# ENCODE long-rna-seq-pipeline (DNAnexus Platform App)

This folder contains the dnanexis applets used in the ENCODE long-rna-seq pipeline. All applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

## Applets (aka steps):
*Preparation:* Typically these steps are only one once per gender specific reference genome (e.g. GCRh37/hg19 female) and
             gene annotation (e.g. gencode.v19).  The result of the preparation steps are used in all subsequent pipeline 
             runs for the same genome, gender and annotation.
- merge-annotation - Takes a gene annotation (e.g. gencode.v19.annotation.gtf.gz) gzipped gtf file, the corresponding tRNA
                     gzipped gtf and a spike-in set (e.g. ERCC) gzipped fasta file and combines them into a single 
                     merged annotation gzipped gtf file.  This file will be used as input to all three 'prep' indexing
                     steps.
- prep-star        - Takes a gender specific genome reference gzipped fasta file (e.g. GCRh37/hg19 female) a merged
                     annotation file and produces a STAR genome index tar.gz file.
- prep-tophat      - Takes a gender specific genome reference gzipped fasta file (e.g. GCRh37/hg19 female) a merged
                     annotation file and produces a TopHat genome index tar.gz file.
- prep-rsem        - Takes a gender specific genome reference gzipped fasta file (e.g. GCRh37/hg19 female) a merged
                     annotation file and produces an RSEM index tar.gz file.

*Normal pipline:* The normal pipeline is actually 2 separate pipelines, as paired-end and single end fastqs are handled
                differently.  The '-pe' and '-se' alignment steps are an obvious distinction, but the pe always produces
                four signal files, while the se pipeline only produces two. 
- align-star-pe        - Takes a pair of (paired-end) gzipped fastq files and a STAR genome index tarred, gzipped file.  
                         This step produces two bams, one aligned to the genome and one aligned to the annotation.
- align-star-se        - Takes a (single-end) gzipped fastq file and the STAR genome index tar.gz file.  
                         This step produces two bams, one aligned to the genome and one aligned to the annotation.
- align-tophat-p       - Takes a pair of (paired-end) gzipped fastq files and a TopHat genome index tarred, gzipped file. 
                         This step produces a single genome aligned bam.
- align-tophat-se      - Takes a (single-end) gzipped fastq file and the TopHat genome index tar.gz file.
                         This step produces a single genome aligned bam.
- bam-to-bw-stranded   - Takes a STAR or TopHat (paired-end) genome-aligned bam file output and a chromosome 
                         name/length file to produce four bigWig 'signal files for easy display in a genome browser.  
                         The alignment signals are filtered by +/- strand and uniquely mapped vs. all mapped reads.
- bam-to-bw-unstranded - Takes a STAR or TopHat (single-end) genome-aligned bam file output and a chromosome 
                         name/length file to produce two bigWig 'signal' files for easy display in a genome browser.  
                         The alignment signals are filtered into uniquely mapped vs. all mapped reads.
- quant-rsem           - Takes a STAR genome aligned bam (single or paired-end) and an RSEM index tarred, gzipped file.
                         This step produces two quantification csv files, one for genes and one for transcripts.

---------
## Flow
*Prepartion pipeline:*
```
INPUTS:  annatation.gtf.gz      genome-gender.fasta.gz   genome-gender.fasta.gz   genome-gender.fasta.gz
         tRNAs.gtf.gz           spike-ins.fasta.gz       spike-ins.fasta.gz       spike-ins.fasta.gz
         spike-ins.fasta.gz     merge-anno.gtf.gz(a)     merge-anno.gtf.gz(a)     merge-anno.gtf.gz(a)
                |                    |                        |                        |
                V                    V                        V                        V
STEPS:   merge-annotation ====> prep-tophat,             prep-star,               prep-rsem
                |                    |                       |
                V                    V                       V
OUTPUTS: merge-anno.gtf.gz(a)   tophat-index.tgz(b)      star-index.tgz(c)        rsem-index.tgz(d)
```

---------
*Paired-end pipeline:*
```
INPUTS:  read1.fq.gz             read1.fq.gz            tophat-pe.bam(e)     star-genome-pe.bam(f)   star-anno-pe.bam(g)
         read2.fq.gz             read2.fq.gz            chrom.sizes          chrom.sizes             rsem-index.tgz(d)
         tophat-index.tar.gz(b)  star-index.tar.gz(c)        |                     |                     |
              |                       |                      |                     |                     |
              V                       V                      V                     V                     V
STEPS:   align-tophat-pe,        align-star-pe   ====>  bam-to-bw-stranded,  bam-to-bw-stranded ===> quant-rsem
              |                       |                      |                     |                     |
              V                       V                      V                     V                     V
OUTPUTS: tophat-pe.bam(e)        star-genome-pe.bam(f)  tophat-uniq-plus.bw   star-uniq-plus.bw      gene-quant.csv
                                 star-anno-pe.bam(g)    tophat-uniq-minus.bw  star-uniq-minus.bw     transcript-quant.csv
                                                        tophat-all-plus.bw    star-all-plus.bw
                                                        tophat-all-minus.bw   star-all-minus.bw
```

---------
*Single-end pipeline:*
```
INPUTS:  reads.fq.gz             reads.fq.gz            tophat-se.bam(h)     star-genome-se.bam(i)      star-anno-se.bam(j)
         tophat-index.tar.gz(b)  star-index.tar.gz(c)   chrom.sizes          chrom.sizes                rsem-index.tgz(d)
              |                       |                      |                     |                        |
              V                       V                      V                     V                        V
STEPS:   align-tophat-se,        align-star-se   ====>  bam-to-bw-unstranded, bam-to-bw-unstranded ===> quant-rsem
              |                       |                      |                     |                        |
              V                       V                      V                     V                        V
OUTPUTS: tophat-se.bam(h)        star-genome-se.bam(i)  tophat-uniq.bw        star-uniq.bw              gene-quant.csv
                                 star-anno-se.bam(j)    tophat-all.bw         star-all.bw               transcript-quant.csv
```


                                                         
