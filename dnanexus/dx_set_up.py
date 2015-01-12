#!/usr/bin/env python
import argparse
import os
import sys
import json

import dxpy
import requests
from dxencode import dxencode as dxencode

ASSAY_TYPE = 'RNA-seq'
ASSAY_TERM_ID = "OBI:0001271"
PROJECT_NAME = 'long-rna-seq-pipeline'


def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Set up DNA Methylation runs on DNA Nexus')

    ap.add_argument('-n', '--number',
                    help='stop after this number of jobs',
                    default='9999999',
                    required=False)

    return ap.parse_args()

def main():
    cmnd = get_args()

    ## resolve projects
    project = dxencode.resolve_project(PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    pid =  project.get_id()

    applet = dxencode.find_applet_by_name('fastqc-exp', pid )
    (AUTHID, AUTHPW, SERVER) = dxencode.processkey('www')
    query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % ASSAY_TERM_ID
    res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
    exps = res.json()['@graph']

    n = 0
    for exp in exps:
        acc = exp['accession']
        if len(exp['replicates']) > 0:
            if exp['replicates'][0]['library'].get('size_range', "") != '>200':
                print "Skipping %s with wrong library size (%s)" % (acc, exp['replicates'][0]['library'].get('size_range', ""))
                #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                continue
            if exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity_units', "") == "cells":
                ncells = float(exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity', 0.0))
                if ncells < 20:
                    print "Skipping %s as single-cell (%s %s)" % (acc, exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity_units', ""), ncells)
                    #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                    continue
            run = applet.run({ "accession": acc}, project=pid)
            print "Running: %s for %s" % (run, acc)
            n = n + 1
            if n > cmnd.number:
                break
        else:
            print "Skipping %s (0 replicates)" % acc

if __name__ == '__main__':
    main()
