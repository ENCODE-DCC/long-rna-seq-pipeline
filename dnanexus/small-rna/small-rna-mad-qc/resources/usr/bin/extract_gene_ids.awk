#! /bin/awk -f 
# Extracts small-rna gene id's from a gencode annotation (gtf) file for use in mad.qc
# usage: extract_gene_ids.awk gencode.v19.annotation.gtf out=srna_gene_ids.txt
#   input files:
#      gencode.v19.annotation.gtf:  Gencode annotation file
#   output file:
#      srna_gene_expr.tsv:   Just the small RNA genes in the same quantification format

BEGIN {OFS="\t"; ix=0; max=0;};
{
    if ($3=="gene" || substr($14,2,length($14)-3)=="tRNAscan") {
    
        g=substr($14,2,length($14)-3); 
        
        if (g=="miRNA" || g=="snoRNA" || g=="snRNA" || g=="tRNAscan") {
            #print substr($10,2,length($10)-3);
            gn=substr($10,2,length($10)-3)
            C[ix++]=gn;
            max=ix 
        }; 
    };
};

END {
    for (ix=0; ix < max; ix++) {
        print C[ix] > out;
    };
};

