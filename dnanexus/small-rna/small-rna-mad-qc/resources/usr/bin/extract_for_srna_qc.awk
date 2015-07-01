#! /bin/awk -f 
# Extracts small-rna gene id's from STAR quantification files for use in mad.qc
# usage: extract_for_srna_qc.awk smallRNA.geneID star_quant.tsv > srna_gene_expr.tsv
## input files:
# smallRNA.geneID:      Small RNA gene IDs extracted from a gencode alignments gtf and saved in the star index archive.
# star_quant.tsv:       STAR generated gene quantification file (tab separated values)
## output file:
# srna_gene_expr.tsv:   Just the small RNA genes in the same quantification format

# example
#
# extract_for_srna_qc.awk smallRNA.geneID ENCSR000AFU_rep1_1_srna_star_quant.tsv > ENCSR000AFU_rep1_1_srna_gene_expr.tsv

    awk 'BEGIN {OFS="\t"} {if (ARGIND==1) {G[$1]=1} else {if (substr($1,1,2)!="N_") {t+=$3} else {print $1,$3; T[$1]=$3}; if ($1 in G) {C[$1                       ]=$3; s+=$3}}} END {for (g in C) {print g,".",".",".",".",".",C[g]/(t+T["N_noFeature"]+T["N_ambiguous"]+T["N_multimapping"])*1000000>"smallRNA.                       expr"}; print "N_notSmallRNA", t-s; print "N_smallRNA",s }'
       $STARgenomeDir/smallRNA.geneID ReadsPerGene.out.tab > smallRNA.info

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
