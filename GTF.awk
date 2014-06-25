#! /bin/awk -f 
# GTF.awk v1.0 creates GTF file from Gencode annotations, Gencode tRNAs, spike-ins fasta
# usage: GTF.awk main.gtf tRNA.gtf spikein.fasta > combined.gtf
## input files:
# main.gtf:      main gtf file, e.g. "gencode.v19.annotation.gtf"
# tRNA.gtf:      gtf file for tRNAs or other annotations - exons only!
# spikein.fasta: fasta file with spike-ins, e.g. "spikes.fixed.fasta"
## output file:
# combined.gtf:  GTF file to be used for genome generation with STAR, RSEM, TopHat


# example
#
# GTF.awk gencode.v19.annotation.gtf gencode.v19.tRNAs.gtf spikes.fixed.fasta > gencode.v19.annotation.tRNA.gtf

BEGIN {FS="\t";OFS="\t"};

{
    if (substr($1,1,1)!="#") {# remove all comments
        if (ARGIND==1) {# first file - print everyting
            print;
        } else if (ARGIND==2) {
            $3="exon"; 
            print;
        } else if (ARGIND==3) {
            if (substr($1,1,1)==">") {
                if (L>0) {
                    print sp,"spikein","exon",1,L,".","+",".","gene_id \"gSpikein_" sp "\"; transcript_id \"tSpikein_" sp "\";";
                    L=0;
                };
                sp=substr($1,2);
            } else {
                L+=length($1);
            };
        };
    };
};


END {
    print sp,"spikein","exon",1,L,".","+",".","gene_id \"gSpikein_" sp "\"; transcript_id \"tSpikein_" sp "\";";
}
