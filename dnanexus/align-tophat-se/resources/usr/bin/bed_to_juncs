#!/usr/bin/env python
# encoding: utf-8
"""
bed_to_juncs.py

Created by Cole Trapnell on 2008-09-19.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
This script converts junctions in BED format produced by TopHat to the 
internal .juncs format for re-use with future runs.

Usage:

    bed_to_juncs.py < junctions.bed
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
            raise Usage(msg)
    
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)

        line_num = 0
        for line in sys.stdin.readlines():
            line = line.strip()
            cols = line.split()
            line_num += 1
            if len(cols) < 12:
                print >> sys.stderr, "Warning: malformed line %d, missing columns" % line_num
                print >> sys.stderr, "\t", line
                continue
            chromosome = cols[0]
            orientation = cols[5]
            block_starts = [int(x) for x in cols[11].split(",")]
            block_sizes = [int(x) for x in cols[10].split(",")]
                
            left_pos = int(cols[1]) + block_starts[0] + block_sizes[0] - 1
            right_pos = int(cols[1]) + block_starts[1] 
            print "%s\t%d\t%d\t%s" % (chromosome, left_pos, right_pos, orientation)
            

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
