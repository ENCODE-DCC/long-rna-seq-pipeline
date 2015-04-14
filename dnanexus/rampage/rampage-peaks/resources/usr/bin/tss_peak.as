table tss_peak
"Rampage peak coverage. Browser extensible data (11 fields)."
    (
    string  chrom;       "Reference sequence chromosome (or contig, scaffold, etc.)"
    uint    chromStart;  "Start position in chromosome"
    uint    chromEnd;    "End position in chromosome"
    string  name;        "Name of item"
    uint    score;       "Score from 0-1000"
    char[1] strand;      "+ or - for strand"
    float   count;       "Count of reads mapping to this peak"
    string  gene_id;     "Gene identifier"
    string  gene_name;   "Gene name" 
    string  tss_id;      "TSS identifier"
    lstring peak_cov;    "base by base read coverage of the peak"
    )
