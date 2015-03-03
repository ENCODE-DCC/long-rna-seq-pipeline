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
    ap = argparse.ArgumentParser(description='Calculate QC numbers for RNA-seq jobs')

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

    applet = dxencode.find_applet_by_name('lrna-qc', pid )
    (AUTHID, AUTHPW, SERVER) = dxencode.processkey('www')
    query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % ASSAY_TERM_ID
    res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
    exps = res.json()['@graph']

    head = ""
    tabfh = open('lrna_qc.tsv','w')
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


                annd = {}
                repd = {}
                fild = {}

                annd = dxfiles.setdefault(assembly,{})
                repd = annd.setdefault(annotation, {})
                f.update({ 'dxid': json.loads(f['notes'])['dx-id'], 'rstr': rstr })
                fild = { f['output_type']: f }
                if repd.has_key(rstr):
                    repd[rstr].update(fild)
                else:
                    repd[rstr] = fild

                if annd.has_key(annotation):
                    annd[annotation].update(repd)
                else:
                    annd[annotation] = repd

                if dxfiles.has_key(assembly):
                    dxfiles[assembly].update(annd)
                else:
                    dxfiles[assembly] = annd

        else:
            print( "Skipping %s (0 replicates)" % (acc))

        r1qs = []
        r2qs = []
        for assembly in dxfiles.keys():
            for annotation in dxfiles[assembly].keys():
                for rep1 in dxfiles[assembly][annotation].keys():
                    for rep2 in [ r for r in dxfiles[assembly][annotation].keys() if r != rep1 ]:
                        for out_type in dxfiles[assembly][annotation][rep1].keys():
                            if out_type.find('quantification') >= 0:
                                r1qs.append(dxfiles[assembly][annotation][rep1][out_type])
                                r2qs.append(dxfiles[assembly][annotation][rep2][out_type])

        if r1qs and r2qs:
            run = "test"
            print("Running: %s for (%s,%s) in %s" % (run, r1qs, r2qs, acc))
            job = applet.run({ "rep1_quants": [ dxpy.dxlink(f['dxid']) for f in r1qs ], "rep2_quants": [ dxpy.dxlink(f['dxid']) for f in  r2qs ]}, project=pid)
            job.wait_on_done(interval=1)
            pair = 0
            for qcs_str in job.describe()['output'].get('qc_metrics_json', []):
                qcs = json.loads(qcs_str)
                if not head:
                    head= "\t".join(['Experiment', 'Assembly', 'Annotation', 'RepA', 'RepB', 'type']+qcs.keys())
                    tabfh.write(head+"\n")

                tabfh.write("\t".join([ r1qs[pair]['dataset'],
                                  r1qs[pair]['assembly'],
                                  r1qs[pair]['genome_annotation'],
                                  r1qs[pair]['rstr'],
                                  r2qs[pair]['rstr'],
                                  r1qs[pair]['output_type'] ]+[ str(v) for v in qcs.values() ])+"\n" ) # they are floats
                pair += 1
                tabfh.flush()
        elif dxfiles:
            print("%s has DXfiles but nothing to run: %s" % (acc, dxfiles))

    tabfh.close()
if __name__ == '__main__':
    main()
