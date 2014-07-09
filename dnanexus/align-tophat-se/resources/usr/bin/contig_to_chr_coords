#!/usr/bin/env python
# encoding: utf-8
"""
contig_to_chr_coords.py

Created by Cole Trapnell on 2008-09-05.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
Takes the NCBI seq_contig file and maps contig coords to whole chromosome 
coords in a GTF, GFF, or BED file

    contig_to_chr_coords.py <format_flag> <seq_contig.md> <junctions.bed|transcripts.gff|transcripts.gtf>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:vbg", ["help", "output=", "bed", "gff"])
        except getopt.error, msg:
            raise Usage(msg)
    
        arg_is_splice = False
        arg_is_gff = False
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            if option in ("-b", "--bed"):
                arg_is_splice = True
            if option in ("-g", "--gff"):
                arg_is_gff = True
        
        if (arg_is_splice == False and arg_is_gff == False) or (arg_is_splice == True and arg_is_gff == True):
			print >> sys.stderr, "Error: please specify either -b or -g, but not both"
			raise Usage(help_message)
        
        if len(args) < 1:
            raise Usage(help_message)
        contig_to_chr_file = open(args[0])
        
        contigs = {}
        
        for line in contig_to_chr_file.readlines():
            if line[0] == "#":
                continue
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            chromosome = cols[1]
            group = cols[8]
            feature_name = cols[5]
            if not feature_name in ["start", "end"]:
                contigs[feature_name] = (chromosome, int(cols[2]))
                #print feature_name, chromosome, int(cols[2])
        
        if arg_is_gff:
            gff_file = open(args[1])
        
            lines = gff_file.readlines()
            print lines[0],
            for line in lines[1:]:
                line = line.strip()
                cols = line.split('\t')
                if len(cols) < 8:
                    continue
                    
                contig = cols[0]
                chr_fields = contig.split('|')
                refseq_id = chr_fields[3]
                contig = contigs.get(refseq_id)
                chr_name = contig[0]
                pipe_idx = chr_name.find('|')

                if pipe_idx != -1:
                    chr_name = chr_name[:pipe_idx]
                if contig != None:
                    #print line
                    left_pos = contig[1] + int(cols[3])
                    right_pos = contig[1] + int(cols[4])
                
                    print "chr%s\tTopHat\tisland\t%d\t%d\t%s\t.\t.\t%s" % (chr_name, left_pos, right_pos, cols[5],cols[8])
                    #print >>sys.stderr, "%s\t%d\t%d\t%s\t%s\t%s\t%s" % (contig[0], left_pos, right_pos,cols[3],cols[6],cols[0],cols[1])
            
        
        if arg_is_splice:
            splice_file = open(args[1])
        
            lines = splice_file.readlines()
            print lines[0],
            for line in lines[1:]:
                line = line.strip()
                cols = line.split('\t')
                contig = cols[0]
                chr_fields = contig.split('|')
                refseq_id = chr_fields[3]
                contig = contigs.get(refseq_id)
                chr_name = contig[0]
                pipe_idx = chr_name.find('|')

                if pipe_idx != -1:
                    chr_name = chr_name[:pipe_idx]
                if contig != None:
                    #print line
                    left_pos = contig[1] + int(cols[1])
                    right_pos = contig[1] + int(cols[2])
                
                    print "chr%s\t%d\t%d\t%s\t0\t%s\t%s\t%s\t255,0,0\t2\t1,1\t%s" % (chr_name, left_pos, right_pos, cols[3],cols[5],left_pos, right_pos,cols[11])
                    #print >>sys.stderr, "%s\t%d\t%d\t%s\t%s\t%s\t%s" % (contig[0], left_pos, right_pos,cols[3],cols[6],cols[0],cols[1])
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
