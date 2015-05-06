#!/usr/bin/env python2.7
# versions.py 0.0.1
#
# Creates versions json string for a particular applet

import sys, os, argparse, json, commands

# APP_TOOLS is a dict keyed by applet script name with a list of tools that it uses.
APP_TOOLS = {
    "prep-rsem.sh":                [ "RSEM" ], 
    "prep-star.sh":                [ "STAR" ], 
    "prep-tophat.sh":              [ "TopHat", "bowtie2" ], 
    "merge-annotation.sh":         [ "GTF.awk" ],
    "concat-fastqs.sh":            [],

    "align-star-pe.sh":            [ "STAR", "samtools" ],
    "align-star-se.sh":            [ "STAR", "samtools" ],
    "align-tophat-pe.sh":          [ "TopHat", "bowtie2", "samtools", "tophat_bam_xsA_tag_fix.pl" ],
    "align-tophat-se.sh":          [ "TopHat", "bowtie2", "samtools" ],
    "bam-to-bigwig-stranded.sh":   [ "STAR","bedGraphToBigWig" ],
    "bam-to-bigwig-unstranded.sh": [ "STAR","bedGraphToBigWig" ],
    "quant-rsem.sh":               [ "RSEM" ],
    
    "small-rna-prep-star.sh":      [ "STAR" ], 
    "small-rna-align.sh":          [ "STAR" ], 
    "small-rna-signals.sh":        [ "STAR","bedGraphToBigWig" ], 

    "rampage-align-pe.sh":         [ "STAR" ],
    "rampage-signals.sh":          [ "STAR", "bedGraphToBigWig" ],
    "rampage-peaks.sh":            [ "call_peaks (grit)", "bedToBigBed", "samtools" ],
    "rampage-idr.sh":              [ "Anaconda3", "idr", "bedToBigBed" ]
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
            "RSEM":                      "rsem-calculate-expression --version | awk '{print $5}'",
            "samtools":                  "samtools 2>&1 | grep Version | awk '{print $2}'",
            "STAR":                      "STAR --version | awk '{print $1}' | cut -d _ -f 2-",
            "TopHat":                    "tophat -v | awk '{print $2}'",
            "tophat_bam_xsA_tag_fix.pl": "perl /usr/bin/tophat_bam_xsA_tag_fix.pl --version 2>&1"
            }

def main():
    parser = argparse.ArgumentParser(description =  "Versions parser for a dx applet. " + \
                                                    "Prints version lines to stderr and json string to stdout.")
    parser.add_argument('-a','--applet', required=True,
                        help="Applet to print versions for")
    parser.add_argument('-av','--appver', required=True,
                        help="Version of applet")
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, 
                        help="Don't print versions to stderr.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Show the command-line that is used to get the version.")


    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 3:
        parser.print_usage()
        return
        
    versions = {}
    versions["DX applet"] = { args.applet: args.appver }
    if not args.quiet:
        sys.stderr.write("********\n")
        sys.stderr.write("* Running " + args.applet + ": " + args.appver+ "\n")
        
    tools = APP_TOOLS[args.applet]
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


