table idr_peak
"IDR peaks for Rampage. Browser extensible data (14 fields)."
    (
    string  chrom;           "Reference sequence chromosome (or contig, scaffold, etc.)"
    uint    chromStart;      "Start position in chromosome"
    uint    chromEnd;        "End position in chromosome"
    string  name;            "Name of item"
    uint    score;           "Score from 0-1000"
    char[1] strand;          "+ or - for strand"
    float   localIDR;        "Local IDR value" 
    float   globalIDR;       "Global IDR value" 
    uint    rep1_chromStart; "Start position in chromosome of replicate 1 peak"    
    uint    rep1_chromEnd;   "End position in chromosome of replicate 1 peak"    
    float   rep1_count;      "Count (used for ranking) replicate 1"
    uint    rep2_chromStart; "Start position in chromosome of replicate 2 peak"    
    uint    rep2_chromEnd;   "End position in chromosome of replicate 2 peak"    
    float   rep2_count;      "Count (used for ranking) replicate 2"
    )
