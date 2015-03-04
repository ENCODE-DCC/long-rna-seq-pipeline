#!/usr/bin/env python
import argparse
import os
import sys
import json
import random

import dxpy
import requests
from dxencode import dxencode as dxencode

ASSAY_TYPE = 'RNA-seq'
ASSAY_TERM_ID = "OBI:0001271"
PROJECT_NAME = 'long-rna-seq-pipeline'
TOTAL_CONTROLS = 10 # this could be a parameter


def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Calculate QC numbers for RNA-seq jobs')

    ap.add_argument('-n', '--number',
                    help='stop after this number of jobs',
                    default='9999999',
                    required=False)

    ap.add_argument('-c', '--controls',
                    help='generate this many negative control pairs',
                    type=int,
                    default=0,
                    required=False)

    return ap.parse_args()

def run_pairs(head, applet, pid, r1qs, r2qs, outfh=sys.stdout):
    job = applet.run({ "rep1_quants": [ dxpy.dxlink(f['dxid']) for f in r1qs ], "rep2_quants": [ dxpy.dxlink(f['dxid']) for f in  r2qs ]}, project=pid)
    job.wait_on_done(interval=1)
    pair = 0
    for qcs_str in job.describe()['output'].get('qc_metrics_json', []):
        qcs = json.loads(qcs_str)
        if not head:
            head= "\t".join(['Experiment', 'Assembly', 'Annotation', 'RepA', 'RepB', 'type']+qcs.keys())
            outfh.write(head+"\n")

        outfh.write("\t".join([ r1qs[pair]['dataset'],
                          r1qs[pair]['assembly'],
                          r1qs[pair]['genome_annotation'],
                          r1qs[pair]['rstr'],
                          r2qs[pair]['rstr'],
                          r1qs[pair]['output_type'] ]+[ str(v) for v in qcs.values() ])+"\n" ) # they are floats
        pair += 1
        outfh.flush()

    return head
    # probably return some object.

def pick(dictionary):

    while(True):
        mypick = random.choice(dictionary.values())
        if mypick:
            return mypick
        else:
            print "Whiff on %s" % dictionary


def main():
    cmnd = get_args()

    ## resolve projects
    project = dxencode.resolve_project(PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    pid =  project.get_id()

    applet = dxencode.find_applet_by_name('lrna-qc', pid )
    (AUTHID, AUTHPW, SERVER) = dxencode.processkey('www')
    query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % ASSAY_TERM_ID
    res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
    exps = res.json()['@graph']

    head = ""
    tabfh = open('lrna_qc.tsv','w')
    byexperiment = {}
    for exp in exps:
        acc = exp['accession']
        if len(exp['replicates']) > 0:
            if exp['replicates'][0]['library'].get('size_range', "") != '>200':
                print "Skipping %s with wrong library size (%s)" % (acc, exp['replicates'][0]['library'].get('size_range', ""))
                #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                continue
            dxfiles = {}
            for facc in exp['original_files']:
                f = dxencode.encoded_get(SERVER+facc, AUTHID=AUTHID, AUTHPW=AUTHPW).json()
                try:
                    json.loads(f['notes'])
                    assembly = f['assembly']
                    annotation = f['genome_annotation']
                    rstr = "%s_%s" % (f['replicate']['biological_replicate_number'], f['replicate']['technical_replicate_number'])
                except Exception, e:
                    print("Skipping %s/%s (%s) 'cause %s" % (acc, f['accession'], f['output_type'], e))
                    continue

                fild = dxfiles.setdefault(assembly, {}).setdefault(annotation, {}).setdefault(rstr, {})
                f.update({ 'dxid': json.loads(f['notes'])['dx-id'], 'rstr': rstr })
                fild[f['output_type']] = f

            byexperiment[acc] = dxfiles

        else:
            print( "Skipping %s (0 replicates)" % (acc))
            continue

        r1qs = []
        r2qs = []

        for assembly in dxfiles.values():
            for annotation in assembly.values():
                for rep1 in annotation.values():
                    for rep2 in annotation.values():
                        if rep2 == rep1:
                            continue
                        for out_type, out_type_value in rep1.items():
                            if out_type.find('quantification') >= 0:
                                r1qs.append(out_type_value)
                                r2qs.append(rep2[out_type])


            if r1qs and r2qs:
                head = run_pairs(head, applet, pid, r1qs, r2qs, outfh=tabfh)
            elif dxfiles:
                print("%s has DXfiles but nothing to run: %s" % (acc, dxfiles))

    tabfh.close()
    if cmnd.controls:
        control1s = []
        control2s = []
        ctfh = open('lrna_qcs_controls.tsv','w')
        head = ''
        for c in range(0, cmnd.controls):
            for qtype in ('transcript quantifications', 'genome quantifications'):
                control1s.append(pick(pick(pick(pick(byexperiment))))[qtype])
                control2s.append(pick(pick(pick(pick(byexperiment))))[qtype])

        run_pairs(head, applet, pid, control1s, control2s, outfh=ctfh)
if __name__ == '__main__':
    main()
