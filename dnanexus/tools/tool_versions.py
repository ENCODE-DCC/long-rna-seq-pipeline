#!/usr/bin/env python2.7
# tool_versions.py v1.1  Creates "SW" versions json string for a particular DX applet.
#                        Write request to stdout and verbose info to stderr.  This allows easy use in dx app scripts.

import sys, os, argparse, json, commands

# APP_TOOLS is a dict keyed by applet script name with a list of tools that it uses.
APP_TOOLS = {
    # lrna:    
    "align-star-pe":            [ "lrna_align_star_pe.sh", "STAR", "samtools" ],
    "align-star-se":            [ "lrna_align_star_se.sh", "STAR", "samtools" ],
    "align-tophat-pe":          [ "lrna_align_tophat_pe.sh", "TopHat", "bowtie2", "samtools", "tophat_bam_xsA_tag_fix.pl" ],
    "align-tophat-se":          [ "lrna_align_tophat_se.sh", "TopHat", "bowtie2", "samtools" ],
    "bam-to-bigwig-stranded":   [ "lrna_bam_to_stranded_signals.sh", "STAR","bedGraphToBigWig" ],
    "bam-to-bigwig-unstranded": [ "lrna_bam_to_unstranded_signals.sh", "STAR","bedGraphToBigWig" ],
    "quant-rsem":               [ "lrna_rsem_quantification.sh", "RSEM" ],
    "mad-qc":                   [ "MAD.R" ],

    # srna:    
    "small-rna-prep-star":      [ "srna_index.sh", "STAR", "extract_gene_ids.awk" ], 
    "small-rna-align":          [ "srna_align.sh", "STAR", "samtools" ], 
    "small-rna-signals":        [ "srna_signals.sh", "STAR","bedGraphToBigWig" ], 
    "small-rna-mad-qc":         [ "srna_mad_qc.sh", "MAD.R", "extract_gene_ids.awk", "sum_srna_expression.awk" ],

    # rampage:    
    "rampage-align-pe":         [ "rampage_align_star.sh", "STAR", "samtools" ],
    "rampage-signals":          [ "rampage_signal.sh", "STAR", "bedGraphToBigWig" ],
    "rampage-peaks":            [ "rampage_peaks.sh", "call_peaks (grit)", "bedToBigBed", "pigz", "samtools" ],
    "rampage-idr":              [ "rampage_idr.sh", "Anaconda3", "idr", "bedToBigBed", "pigz" ],
    "rampage-mad-qc":           [ "rampage_mad_qc.sh", "MAD.R" ],

    # utility:    
    "merge-annotation":         [ "GTF.awk" ],
    "prep-rsem":                [ "lrna_index_rsem.sh", "RSEM" ], 
    "prep-star":                [ "lrna_index_star.sh", "STAR" ], 
    "prep-tophat":              [ "lrna_index_tophat.sh", "TopHat", "bowtie2" ], 
    }
# Virtual apps only differ from their parent by name/version. 
VIRTUAL_APPS = {
    # lrna virtuals:    
    "bam-to-bigwig-stranded-tophat":   "bam-to-bigwig-stranded",
    "bam-to-bigwig-unstranded-tophat": "bam-to-bigwig-unstranded",
    "quant-rsem-alt":                  "quant-rsem",
    "mad-qc-alt":                      "mad-qc",

    }


# ALL_TOOLS contains the printable tool name (key) and the command that is used to determine the version.
ALL_TOOLS = {
            "Anaconda3":                 "ls Anaconda3*.sh | head -1 | cut -d - -f 2",
            "bedGraphToBigWig":          "bedGraphToBigWig 2>&1 | grep 'bedGraphToBigWig v' | awk '{print $2$3}'",
            "bedToBigBed":               "bedToBigBed 2>&1 | grep 'bedToBigBed v' | awk '{print $2$3}'",
            "bowtie2":                   "bowtie2 --version 2>&1 | grep bowtie | awk '{print $3}'",
            "call_peaks (grit)":         "call_peaks --version 2>&1 | grep call_peaks | awk '{print $3}'",
            "GTF.awk":                   "echo unversioned",
            "idr":                       "idr/bin/idr --version 2>&1 | grep IDR | awk '{print $2}'",
            "MAD.R":                     "grep version /usr/bin/MAD.R | awk '{print $3}'",
            "extract_gene_ids.awk":      "grep version /usr/bin/extract_gene_ids.awk | awk '{print $3}'",
            "sum_srna_expression.awk":   "grep version /usr/bin/sum_srna_expression.awk | awk '{print $3}'",
            "RSEM":                      "rsem-calculate-expression --version | awk '{print $5}'",
            "samtools":                  "samtools 2>&1 | grep Version | awk '{print $2}'",
            "STAR":                      "STAR --version | awk '{print $1}' | cut -d _ -f 2-",
            "TopHat":                    "tophat -v | awk '{print $2}'",
            "tophat_bam_xsA_tag_fix.pl": "perl /usr/bin/tophat_bam_xsA_tag_fix.pl --version 2>&1",
            "pigz":                      "pigz --version 2>&1 | awk '{print $2}'",
            "lrna_align_star_pe.sh":             "lrna_align_star_pe.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_align_star_se.sh":             "lrna_align_star_se.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_align_tophat_pe.sh":           "lrna_align_tophat_pe.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_align_tophat_se.sh":           "lrna_align_tophat_se.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_bam_to_stranded_signals.sh":   "lrna_bam_to_stranded_signals.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_bam_to_unstranded_signals.sh": "lrna_bam_to_unstranded_signals.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_rsem_quantification.sh":       "lrna_rsem_quantification.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_index_rsem.sh":                "lrna_index_rsem.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_index_star.sh":                "lrna_index_star.sh | grep usage | awk '{print $2}' | tr -d :",
            "lrna_index_tophat.sh":              "lrna_index_tophat.sh | grep usage | awk '{print $2}' | tr -d :",
            "rampage_align_star.sh":    "rampage_align_star.sh | grep usage | awk '{print $2}' | tr -d :",
            "rampage_signal.sh":        "rampage_signal.sh | grep usage | awk '{print $2}' | tr -d :",
            "rampage_peaks.sh":         "rampage_peaks.sh | grep usage | awk '{print $2}' | tr -d :",
            "rampage_idr.sh":           "rampage_idr.sh | grep usage | awk '{print $2}' | tr -d :",
            "rampage_mad_qc.sh":        "rampage_mad_qc.sh | grep usage | awk '{print $2}' | tr -d :",
            "srna_index.sh":            "srna_index.sh | grep usage | awk '{print $2}' | tr -d :",
            "srna_align.sh":            "srna_align.sh | grep usage | awk '{print $2}' | tr -d :",
            "srna_signals.sh":          "srna_signals.sh | grep usage | awk '{print $2}' | tr -d :",
            "srna_mad_qc.sh":           "srna_mad_qc.sh | grep usage | awk '{print $2}' | tr -d :",
            }

def parse_dxjson(dxjson):
    '''Parses the dnanexus-executable.json file in the job directory to get applet name and version.'''
    with open(dxjson) as data_file:    
        dxapp = json.load(data_file)

    appver = "unknown"    
    applet = dxapp.get("name").split()[0]
    if "version" in dxapp:
        appver = dxapp.get("version")
    else:
        title = dxapp.get("title")
        last_word = title.split(' ')[-1]
        if last_word.startswith('(virtual-') and last_word.endswith(')'):
            appver = last_word[9:-1]
        elif last_word.startswith('(v') and last_word.endswith(')'):
            appver = last_word[2:-1]
    
    return (applet, appver)


def main():
    parser = argparse.ArgumentParser(description =  "Versions parser for a dx applet. " + \
                                                    "Prints version lines to stderr and json string to stdout. " + \
                                                    "MUST specify either --applet and --appver or --dxjson.")
    parser.add_argument('-a','--applet', required=False,
                        help="Applet to print versions for")
    parser.add_argument('-av','--appver', required=False,
                        help="Version of applet")
    parser.add_argument('-j','--dxjson', required=False,
                        help="Use dnanexus json file to discover 'applet' and 'appver'")
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, 
                        help="Don't print versions to stderr.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Show the command-line that is used to get the version.")


    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 3:
        parser.print_usage()
        return
        
    if (args.applet == None or args.appver == None) and args.dxjson == None:
        parser.print_help()
        return

    applet = args.applet
    applet = args.appver
    
    if args.dxjson != None:
        (applet,appver) = parse_dxjson(args.dxjson)
    
    versions = {}
    versions["DX applet"] = { applet: appver }
    if not args.quiet:
        sys.stderr.write("********\n")
        sys.stderr.write("* Running " + applet + ": " + appver+ "\n")
    
    if applet in VIRTUAL_APPS:
        tools = APP_TOOLS[VIRTUAL_APPS[applet]]
    else:
        tools = APP_TOOLS[applet]
    for tool in tools:
        cmd = ALL_TOOLS[tool]
        if args.verbose:
            sys.stderr.write("cmd> " + cmd + "\n")
        err, ver = commands.getstatusoutput(cmd)
        versions[tool] = ver
        if not args.quiet:
            sys.stderr.write("* " + tool + " version: " + ver + "\n")

    if not args.quiet:
        sys.stderr.write("********\n")
    
    print json.dumps(versions) 
     
if __name__ == '__main__':
    main()


