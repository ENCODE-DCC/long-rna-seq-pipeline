#!/usr/bin/env python
import argparse
import os

from dxencode import dxencode as dxencode

ASSAY_TYPE = 'RNA-seq'
ASSAY_TERM_ID = "OBI:0001271"
PROJECT_NAME = 'long-rna-seq-pipeline'

GENOME_MAPPING = {
    "human": "hg19",
    "mouse": "mm10"

}

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Set up DNA Methylation runs on DNA Nexus')

    ap.add_argument('-t', '--test',
                    help='Use test input folder',
                    action='store_true',
                    required=False)

    ap.add_argument('-n', '--numberjobs',
                    help='stop after this number of jobs',
                    default=9999999,
                    type=int,
                    required=False)


    return ap.parse_args()

def main():
    cmnd = get_args()

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
        if n >= cmnd.numberjobs:
            print "Stopping at %s replicates" % n
            break
        exp_mapping = dxencode.choose_mapping_for_experiment(exp)
        for rep in exp.get('replicates', []):
            try:
                br = rep['biological_replicate_number']
                tr = rep['technical_replicate_number']
                mapping = exp_mapping[(br,tr)]
                o = GENOME_MAPPING[mapping['organism']]
                args = "-o %s" % o
                args += " -l %s" % mapping['library']
                args += " -g %s" % mapping['sex']
                if mapping['paired']:
                    paired_fqs = {
                        '1': [],
                        '2': []
                    }
                    for (p1, p2) in mapping['paired']:
                        paired_fqs[p1['paired_end']].append(p1['accession']+".fastq.gz")
                        paired_fqs[p2['paired_end']].append(p2['accession']+".fastq.gz")
                    args += " -1 " + " ".join(paired_fqs['1'])
                    args += " -2 " + " ".join(paired_fqs['2'])
                else:
                    args += " -1 " + " ".join([ f['accession']+".fastq.gz" for f in mapping['unpaired'] ])

                runcmd = "./lrnaLaunch.py -e %s -r %s -tr %s %s -a M4 --project %s --resultsLoc /runs --run > runs/launch%s-%s-%s-M4.%s.out" % (acc, br, tr, args, PROJECT_NAME, acc, br, tr, os.getpid())
                print runcmd
                if not cmnd.test:
                    # probably should be subprocess.Popen()
                    os.system(runcmd)
                n+=1
            except KeyError, e:
                print "%s failed: %s" % (acc, e)

if __name__ == '__main__':
    main()
