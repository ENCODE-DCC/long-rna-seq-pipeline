#! /bin/awk -f 
# Sums small-rna gene expression from STAR quantification files for use in mad.qc
# version 1.0
# usage: sum_srna_expression.awk smallRNA.geneID star_quant.tsv out=srna_gene_expr.tsv > srna_expr_summary.txt
#    input files:
#       smallRNA.geneID:        Small RNA gene IDs extracted from a gencode alignments gtf and saved in the star index archive.
#       star_quant.tsv:         STAR generated gene quantification file (tab separated values)
#    output file:
#       srna_gene_expr.tsv:     Small RNA gene expression
#       srna_expr_summary.txt:  Summary of expression

BEGIN {OFS="\t"};

{
    if (ARGIND==1) {
        G[$1]=1;
    } else {
        if (substr($1,1,2)!="N_") {
            t+=$3;
        } else {
            print $1,$3; 
            T[$1]=$3;
        }; 
        if ($1 in G) {
            C[$1]=$3; 
            s+=$3;
        };
    };
};

END {
    for (g in C) {
        print g,".",".",".",".",".",C[g]/(t+T["N_noFeature"]+T["N_ambiguous"]+T["N_multimapping"])*1000000>out;
    }; 
    print "N_notSmallRNA", t-s; 
    print "N_smallRNA",s; 
};

