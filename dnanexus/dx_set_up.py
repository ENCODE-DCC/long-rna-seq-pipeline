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
    query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.size_range=>200&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % ASSAY_TERM_ID
    res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
    exps = res.json()['@graph']

    import pdb;pdb.set_trace()
    n = 0
    for exp in exps:
        acc = exp['accession']
        if len(exp['replicates']) > 0:
            if exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity_units', "") == "cells":
                print "Skipping %s as single-cell" % acc
                print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
            #run = applet.run({ "accession": acc}, project=pid)
            run = "test"
            print "Running: %s for %s" % (run, acc)
            n = n + 1
            if n > cmnd.number:
                break
        else:
            print "Skipping %s (0 replicates)" % acc

if __name__ == '__main__':
    main()
