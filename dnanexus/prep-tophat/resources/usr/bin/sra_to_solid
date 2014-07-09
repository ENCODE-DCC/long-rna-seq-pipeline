#!/usr/bin/env python

"""
sra_to_solid.py
"""

import sys

use_message = '''
convert SOLiD sequences downloaded from SRA FTP (not via the web interface) to a format TopHat processes.
the script simply removes one primer quality value '!' from the sequences.

Usage:
    sra_to_solid input.fastq > output.fastq
'''

if __name__ == "__main__":
    if len(sys.argv) == 2:
        input_file = open(sys.argv[-1], 'r')
        expect_qual = 0
        for line in input_file:
            line = line.rstrip('\n')
            if expect_qual % 4  == 3:
                line = line[1:]
            
            print line
            expect_qual = (expect_qual + 1) % 4
                
    else:
        print use_message;
