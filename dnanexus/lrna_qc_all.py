#!/usr/bin/env python
import argparse
import os
import sys
import json
import random

import dxpy
import dxpy.exceptions
import requests
from dxencode import dxencode as dxencode

ASSAY_TERM_NAME = 'RNA-seq'
ASSAY_TERM_ID = "OBI:0001271"
PROJECT_NAME = 'long-rna-seq-pipeline'


def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Calculate QC numbers for RNA-seq jobs')

    ap.add_argument('-c', '--controls',
                    help='generate this many negative control pairs',
                    type=int,
                    default=0,
                    required=False)

    ap.add_argument('-o','--only-controls',
                    help='Skip the real qc checks, just do controls',
                    action='store_true',
                    required=False)

    ap.add_argument('-r,','--reps',
                    help='Replicates of the form ENCSR###,biorep,techrep ',
                    nargs='+',
                    required=False)

    ap.add_argument('-g,', '--assembly',
                    help='Genome assembly of the selected pair',
                    default='mm10',
                    required=False)

    ap.add_argument('-a,', '--annotation',
                    help='Transcript (Annotation) of the selected replicate pair',
                    default='M4',
                    required=False)

    ap.add_argument('-t,', '--type',
                    help='Use only genome or transcript files',
                    choices=['g','t','genome','transcript','genome quantifications', 'transcript quantification'],
                    default='genome quantifications',
                    required=False)

    return ap.parse_args()

def run_pairs(head, applet, pid, r1qs, r2qs, outfh=sys.stdout):
    try:
        job = applet.run({ "rep1_quants": [ dxpy.dxlink(f['dxid']) for f in r1qs ], "rep2_quants": [ dxpy.dxlink(f['dxid']) for f in  r2qs ]}, project=pid)
        job.wait_on_done(interval=1)
    except dxpy.exceptions.DXJobFailureError, e:
        print("ERROR: %s" % (e))
        return ""

    pair = 0
    for qcs_str in job.describe()['output'].get('qc_metrics_json', []):
        qcs = json.loads(qcs_str)
        if not head:
            head= "\t".join(['Experiment', 'Assembly', 'Annotation', 'RepA', 'RepB', 'type']+qcs.keys())
            outfh.write(head+"\n")

        labels = {}
        for field in ('dataset', 'assembly','genome_annotation','output_type'):
            if r1qs[pair][field] != r2qs[pair][field]:
                labels[field] = '/'.join((r1qs[pair][field].strip('/experiments/'), r2qs[pair][field].strip('/experiments/')))
            else:
                labels[field] = r1qs[pair][field]

        outfh.write("\t".join([ labels['dataset'],
                          labels['assembly'],
                          labels['genome_annotation'],
                          r1qs[pair]['rstr'],
                          r2qs[pair]['rstr'],
                          labels['output_type'] ]+[ str(v) for v in qcs.values() ])+"\n" ) # they are floats
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

def get_file(acc, facc):
    f = dxencode.encoded_get(SERVER+facc, AUTHID=AUTHID, AUTHPW=AUTHPW).json()
    try:
        json.loads(f['notes'])
        assembly = f['assembly']
        annotation = f['genome_annotation']
        rstr = "%s_%s" % (f['replicate']['biological_replicate_number'], f['replicate']['technical_replicate_number'])
    except Exception, e:
        print("Skipping %s/%s (%s) 'cause %s" % (acc, f['accession'], f['output_type'], e))
        return({},'','','')

    f.update({ 'dxid': json.loads(f['notes'])['dx-id'], 'rstr': rstr })

    return f,rstr,assembly, annotation


(AUTHID, AUTHPW, SERVER) = dxencode.processkey('www')

def main():
    cmnd = get_args()

    ## resolve projects
    project = dxencode.resolve_project(PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    pid =  project.get_id()

    applet = dxencode.find_applet_by_name('lrna-qc', pid )
    head = ""

    byexperiment = {}
    r1qs = []
    r2qs = []

    if cmnd.reps:
        if len(cmnd.reps) != 2:
            print("ERROR: Please select exactly 2 replicates with -r/--reps or leave blankf for all")
            sys.exit(1)
        n = 1
        tabfh = open('lrna_qc_%s.tsv' % (cmnd.reps[0].replace(',','_')+cmnd.reps[1].replace(',','_')) ,'w')

        for repstr in cmnd.reps:
            (exp_acc,br,tr) = repstr.split(',')
            print("Looking for %s %s %s (%s %s)" % (exp_acc,"%s_%s" %(br,tr), cmnd.type, cmnd.assembly, cmnd.annotation))
            query = '/experiments/%s' % exp_acc
            res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
            exp = res.json()
            if exp['assay_term_id'] != ASSAY_TERM_ID:
                print("ERROR: %s is not an %s experiment" % (exp_acc, ASSAY_TERM_NAME))
                sys.exit(1)
            if not len(exp['replicates']):
                print("ERROR: %s has no replicates" % (exp_acc))

            for facc in exp['original_files']:
                (f, rstr, assembly, annotation) = get_file(exp_acc, facc)
                if f and assembly == cmnd.assembly and annotation == cmnd.annotation and rstr == "%s_%s" % (br,tr):
                    if f['output_type'].find(cmnd.type) >= 0:
                        if n==1:
                            r1qs.append(f)
                        elif n==2:
                            r2qs.append(f)
                        else:
                            print("Something went wrong")
                            sys.exit(1)
                        n+=1
                        print("WARN: Found %s %s %s (%s %s)" % (facc, rstr, f['output_type'], assembly, annotation))
                        break
                elif f:
                    print("WARN: Rejecting %s %s %s (%s %s)" % (facc, rstr, f['output_type'], assembly, annotation))


    else:
        tabfh = open('lrna_qc_all.tsv','w')

        query = 'search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % ASSAY_TERM_ID
        res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
        exps = res.json()['@graph']


        for exp in exps:
            acc = exp['accession']
            if len(exp['replicates']) > 0:
                if exp['replicates'][0]['library'].get('size_range', "") != '>200':
                    print("WARN: Skipping %s with wrong library size (%s)" % (acc, exp['replicates'][0]['library'].get('size_range', "")))
                    #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                    continue
                dxfiles = {}
                for facc in exp['original_files']:
                    (f,rstr,assembly,annotation) = get_file(acc,facc)
                    if not f:
                        continue

                    fild = dxfiles.setdefault(assembly, {}).setdefault(annotation, {}).setdefault(rstr, {})
                    fild[f['output_type']] = f

                byexperiment[acc] = dxfiles

            else:
                print( "Skipping %s (0 replicates)" % (acc))
                continue

            if not cmnd.only_controls:
                for assembly in dxfiles.values():
                    for annotation in assembly.values():
                        for rep1 in annotation.values():
                            for rep2 in annotation.values():
                                if rep1 is rep2:
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
