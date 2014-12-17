#!/usr/bin/env python
import argparse
import os
import sys
import json
import itertools

import dxpy
#import dxencode as dxencode
from dxencode import dxencode as dxencode

# NOTES: This command-line utility will run the long RNA-seq pipeline for a single replicate
#      - All results will be written to a folder /lrna/<expId>/rep<#>.
#      - If any results are already in that directory, then the steps that created those results
#        will not be rerun.
#      - If any jobs for the experiment and replicate are already running, nothing new will be
#        launched.
#      - Most of the code is generic and found in dxencode.  It relies upon hard-coded JSON below.
#        Tokens are used to abstract dx.app input/outout file names to avoid collisions.
#        - STEP_ORDER is the list of steps in the pipeline
#        - STEPS contains step definitions and enforces dependencies by input file tokens matching
#          to result file tokens of earlier steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.

GENOME_DEFAULT = 'hg19'
''' This the default Genome that long RNA-seq experiments are mapped to.'''

ANNO_DEFAULT = 'v19'

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''

RESULT_FOLDER_DEFAULT = '/lrna'
''' This the default location to place results folders for each experiment.'''

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

STEP_ORDER = {
    # for SE or PE the list in order of steps to run
    "se": [ "concatR1",             "align-tophat-se", "topBwSe", "align-star-se", "starBwSe", "quant-rsem" ],
    "pe": [ "concatR1", "concatR2", "align-tophat-pe", "topBwPe", "align-star-pe", "starBwPe", "quant-rsem" ]
    # examples for testing build_simple_steps
    #"se": [ "align-tophat-se", "align-star-se", "quant-rsem" ],
    #"pe": [ "align-tophat-pe", "align-star-pe", "quant-rsem" ]
    }
'''The (artifically) linear order of all pipeline steps for single or paired-end.'''

STEPS = {
    # for each step: app, list of any params, inputs and results (both as fileToken: app_obj_name)
    # TODO: Any results files not listed here would not be 'deprecated' on reruns.
    "concatR1": {
                "app":     "concat-fastqs",
                "params":  { "rootR1": "outfile_root" },
                "inputs":  { "reads1_set": "reads_set" },
                "results": { "reads1": "reads" }
    },
    "concatR2": {
                "app":     "concat-fastqs",
                "params":  { "rootR2": "outfile_root" },
                "inputs":  { "reads2_set": "reads_set" },
                "results": { "reads2": "reads" }
    },
    "align-tophat-se": {
                "app":     "align-tophat-se",
                "params":  { "library_id":   "library_id" }, #, "nthreads"
                "inputs":  { "reads1":       "reads",
                             "tophat_index": "tophat_index" },
                "results": { "tophat_bam":   "tophat_bam" }
    },
    "align-tophat-pe": {
                "app":     "align-tophat-pe",
                "params":  { "library_id":   "library_id" }, #, "nthreads"
                "inputs":  { "reads1":       "reads_1",
                             "reads2":       "reads_2",
                             "tophat_index": "tophat_index" },
                "results": { "tophat_bam":   "tophat_bam" }
    },
    "topBwSe":  {
                "app":     "bam-to-bigwig-unstranded",
                "inputs":  { "tophat_bam":     "bam_file",
                             "chrom_sizes":    "chrom_sizes" },
                "results": { "tophat_all_bw":  "all_bw",
                             "tophat_uniq_bw": "uniq_bw" }
    },
    "topBwPe":  {
                "app":     "bam-to-bigwig-stranded",
                "inputs":  { "tophat_bam":           "bam_file",
                             "chrom_sizes":          "chrom_sizes" },
                "results": { "tophat_minus_all_bw":  "minus_all_bw",
                             "tophat_minus_uniq_bw": "minus_uniq_bw",
                             "tophat_plus_all_bw":   "plus_all_bw",
                             "tophat_plus_uniq_bw":  "plus_uniq_bw" }
    },
    "align-star-se":   {
                "app":     "align-star-se",
                "params":  { "library_id":      "library_id" }, #, "nthreads"
                "inputs":  { "reads1":          "reads",
                             "star_index":      "star_index" },
                "results": { "star_genome_bam": "star_genome_bam",
                             "star_anno_bam":   "star_anno_bam",
                             "star_log":        "star_log" }
    },
    "align-star-pe":   {
                "app":     "align-star-pe",
                "params":  { "library_id":      "library_id" }, #, "nthreads"
                "inputs":  { "reads1":          "reads_1",
                             "reads2":          "reads_2",
                             "star_index":      "star_index" },
                "results": { "star_genome_bam": "star_genome_bam",
                             "star_anno_bam":   "star_anno_bam",
                             "star_log":        "star_log" }
    },
    "starBwSe": {
                "app":     "bam-to-bigwig-unstranded",
                "inputs":  { "star_genome_bam": "bam_file",
                             "chrom_sizes":     "chrom_sizes" },
                "results": { "star_all_bw":     "all_bw",
                             "star_uniq_bw":    "uniq_bw" }
    },
    "starBwPe": {
                "app":     "bam-to-bigwig-stranded",
                "inputs":  { "star_genome_bam":    "bam_file",
                             "chrom_sizes":        "chrom_sizes" },
                "results": { "star_minus_all_bw":  "minus_all_bw",
                             "star_minus_uniq_bw": "minus_uniq_bw",
                             "star_plus_all_bw":   "plus_all_bw",
                             "star_plus_uniq_bw":  "plus_uniq_bw" }
    },
    "quant-rsem":     {
                "app":     "quant-rsem",
                "params":  { "paired":            "paired" },  #, "nthreads", "rnd_seed"
                "inputs":  { "star_anno_bam":     "star_anno_bam",
                             "rsem_index":        "rsem_index" },
                "results": { "rsem_iso_results":  "rsem_iso_results",
                             "rsem_gene_results": "rsem_gene_results" }
                }
    }

FILE_GLOBS = {
    # For looking up previous result files, use wild-cards
    "reads1":               "/*_concatR1.fq.gz",
    "reads2":               "/*_concatR2.fq.gz",
    "tophat_bam":           "/*_tophat.bam",
    "tophat_minus_all_bw":  "/*_tophat_minusAll.bw",
    "tophat_minus_uniq_bw": "/*_tophat_minusUniq.bw",
    "tophat_plus_all_bw":   "/*_tophat_plusAll.bw",
    "tophat_plus_uniq_bw":  "/*_tophat_plusUniq.bw",
    "tophat_all_bw":        "/*_tophat_all.bw",
    "tophat_uniq_bw":       "/*_tophat_uniq.bw",
    "star_genome_bam":      "/*_star_genome.bam",
    "star_anno_bam":        "/*_star_anno.bam",
    "star_log":             "/*_Log.final.out",
    "star_minus_all_bw":    "/*_star_genome_minusAll.bw",
    "star_minus_uniq_bw":   "/*_star_genome_minusUniq.bw",
    "star_plus_all_bw":     "/*_star_genome_plusAll.bw",
    "star_plus_uniq_bw":    "/*_star_genome_plusUniq.bw",
    "star_all_bw":          "/*_star_genome_all.bw",
    "star_uniq_bw":         "/*_star_genome_uniq.bw",
    "rsem_iso_results":     "/*_rsem.isoforms.results",
    "rsem_gene_results":    "/*_rsem.genes.results"
    }


POST_TEMPLATES = {
    # For looking up previous result files, use wild-cards
    "tophat_bam":   {
        "file_format": "bam",
        "output_type": "alignments",
        "derived_from": ["reads1", "reads2"]
    },
    "tophat_minus_all_bw":  {
        "file_format": "bigWig",
        "output_type": "multi-read minus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_minus_uniq_bw":{
        "file_format": "bigWig",
        "output_type": "unique read minus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_plus_all_bw":   {
        "file_format": "bigWig",
        "output_type": "multi-read plus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_plus_uniq_bw":  {
        "file_format": "bigWig",
        "output_type": "unique read minus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_all_bw":        {
        "file_format": "bigWig",
        "output_type": "multi-read signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_uniq_bw":       {
        "file_format": "bigWig",
        "output_type": "unique signal",
        "derived_from": ["tophat_bam"]
    },
    "star_genome_bam":      {
        "file_format": "bam",
        "output_type": "alignments",
        "derived_from": ["reads1", "reads2"]
    },
    "star_anno_bam":        {
        "file_format": "bam",
        "output_type": "transcriptome alignments",
        "derived_from": ["reads1", "reads2"]
    },
    "star_minus_all_bw":    {
        "file_format": "bigWig",
        "output_type": "multi-read minus signal",
        "derived_from": ["star_genome_bam"]
    },
    "star_minus_uniq_bw":   {
        "file_format": "bigWig",
        "output_type": "unique minus signal",
        "derived_from": ["star_genome_bam"]
    },
    "star_plus_all_bw":     {
        "file_format": "bigWig",
        "output_type": "multi-read plus signal",
        "derived_from": ["star_genome_bam"]
    },
    "star_plus_uniq_bw":    {
        "file_format": "bigWig",
        "output_type": "unique plus signal",
        "derived_from": ["star_genome_bam"]
    },
    "star_all_bw":          {
        "file_format": "bigWig",
        "output_type": "multi-read signal",
        "derived_from": ["star_genome_bam"]
    },
    "star_uniq_bw":         {
        "file_format": "bigWig",
        "output_type": "unique signal",
        "derived_from": ["star_genome_bam"]
    },
    "rsem_iso_results":     {
        "file_format": "tsv",
        "output_type": "transcript quantifications",
        "derived_from": ["star_anno_bam"]
    },
    "rsem_gene_results":    {
        "file_format": "tsv",
        "output_type": "genome quantifications",
        "derived_from": ["star_anno_bam"]
    }
}

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should remove annotation if only one per genome
    # TODO: should use ACCESSION based fileNames
    "tophat_index":  {
                    "hg19": {
                            "female":   {
                                        "v19": "hg19_female_v19_ERCC_tophatIndex.tgz"
                                        },
                            "male":     {
                                        "v19": "hg19_male_v19_ERCC_tophatIndex.tgz"
                                        }
                            },
                    "mm10": {
                            "female":   {
                                        "M2":  "mm10_female_M2_ERCC_tophatIndex.tgz",
                                        "M3":  "mm10_female_M3_ERCC_tophatIndex.tgz",
                                        "M4":  "mm10_female_M4_ERCC_tophatIndex.tgz"
                                        },
                            "male":     {
                                        "M2":  "mm10_male_M2_ERCC_tophatIndex.tgz",
                                        "M3":  "mm10_male_M3_ERCC_tophatIndex.tgz",
                                        "M4":  "mm10_male_M4_ERCC_tophatIndex.tgz"
                                        }
                            }
                    },
    "star_index":    {
                    "hg19": {
                            "female":   {
                                        "v19": "hg19_female_v19_ERCC_starIndex.tgz"
                                        },
                            "male":     {
                                        "v19": "hg19_male_v19_ERCC_starIndex.tgz"
                                        }
                            },
                    "mm10": {
                            "female":   {
                                        "M2":  "mm10_female_M2_ERCC_starIndex.tgz",
                                        "M3":  "mm10_female_M3_ERCC_starIndex.tgz",
                                        "M4":  "mm10_female_M4_ERCC_starIndex.tgz"
                                        },
                            "male":     {
                                        "M2":  "mm10_male_M2_ERCC_starIndex.tgz",
                                        "M3":  "mm10_male_M3_ERCC_starIndex.tgz",
                                        "M4":  "mm10_male_M4_ERCC_starIndex.tgz"
                                        }
                            }
                    },
    "rsem_index":    {
                    "hg19": {
                            "v19": "hg19_male_v19_ERCC_rsemIndex.tgz"
                            },
                    "mm10": {
                            "M2":  "mm10_male_M2_ERCC_rsemIndex.tgz",
                            "M3":  "mm10_male_M3_ERCC_rsemIndex.tgz",
                            "M4":  "mm10_male_M4_ERCC_rsemIndex.tgz"
                            }
                    },
    "chrom_sizes":   {
                    "hg19": {
                            "female":   "female.hg19.chrom.sizes",
                            "male":     "male.hg19.chrom.sizes"
                            },
                    "mm10": {
                            "female":   "female.mm10.chrom.sizes",
                            "male":     "male.mm10.chrom.sizes"
                            }
                    }
    }

GENOME_MAPPING = {
    "human": "hg19",
    "mouse": "mm10"

}
APPLETS = {}
# Used for caching applets that might be called more than once in pipeline
FILES = {}
# Used for caching file dxlinks that might be needed more than once in building the workflow

def get_args():
    '''Parse the input arguments.'''
    ### PIPELINE SPECIFIC
    ap = argparse.ArgumentParser(description="Launches long RNA-seq pipeline analysis for " +
                "one replicate on single or paired-end reads. Can be run repeatedly and will " +
                "launch only the steps that are needed to finish the pipeline. All results " +
                "will be placed in the folder /<resultsLoc>/<experiment>/<replicate>.")
    ### PIPELINE SPECIFIC

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    required=True)

    ap.add_argument('--br', '--biorep',
                    help="Biological Replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=False)

    ap.add_argument('--tr', '--techrep',
                    help="Technical replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=False)


    ### PIPELINE SPECIFIC
    # TODO: should remove annotation if only one per genome
    ap.add_argument('-a', '--annotation',
                    help="Label of annotation (default: '" + ANNO_DEFAULT + "')",
                    choices=[ANNO_DEFAULT, 'M2','M3','M4'],
                    default=ANNO_DEFAULT,
                    required=False)
    ### PIPELINE SPECIFIC

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    default=PROJECT_DEFAULT,
                    required=False)

    ap.add_argument('--refLoc',
                    help="The location to find reference files (default: '" + \
                                            REF_PROJECT_DEFAULT + ":" + REF_FOLDER_DEFAULT + "')",
                    default=REF_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--test',
                    help='Test run only, do not launch anything.',
                    action='store_true',
                    required=False)
    return ap.parse_args()

def pipeline_specific_vars(args, mapping, pairedEnd):
    '''Adds pipeline specific variables to a dict, for use building the workflow.'''
    # psv can contain any variables, but it must contain these at a minimum:
    # - Any non-file input param needed to launch the workflow
    # - 'resultFolder' - full dx path (without project) to the results folder of the specific run
    # - 'name' - A short name used for specific workflow run.
    # - 'description' - generic description of the pipeline.
    # - 'title'/['subtitle'] for command line output announcing what will be done
    # - Should also contain such things as:
    # - 'organism', 'gender', 'experiment', 'replicate' (if appropriate),
    # - 'pairedEnd' (boolean, if appropriate)

    psv = {}
    psv['project']    = args.project
    psv['organism']   = mapping['genome']
    psv['gender']     = mapping['sex']
    psv['annotation'] = args.annotation
    if psv['organism'] == 'hg19' and psv['annotation'] not in [ANNO_DEFAULT]:
        print psv['organism']+" has no "+psv['annotation']+" annotation."
        sys.exit(1)
    if psv['organism'] == 'mm10' and psv['annotation'] not in ['M2','M3','M4']:
        print psv['organism']+" has no '"+psv['annotation']+"' annotation."
        sys.exit(1)
    psv['experiment'] = args.experiment
    psv['replicate']  = str(args.br)
    psv['rep_tech']   = mapping['replicate']
    psv['library_id'] = mapping['library']
    psv['nthreads']   = 8
    psv['rnd_seed']   = 12345
    psv['paired']  = pairedEnd

    ## below is not necessary, strictly
    # workflow labeling
    psv['description'] = "The ENCODE RNA Seq pipeline for long RNAs"
    psv['name'] = "lrna_"+psv['organism']
    if psv['organism'] == 'mm10':
        psv['name'] += psv['annotation']
    if psv['gender'] == 'female':
        psv['name'] += "XX"
    else:
        psv['name'] += "XY"
    if pairedEnd:
        psv['title'] = "long RNA-seq paired-end "
        psv['name'] += "PE"
    else:
        psv['title'] = "long RNA-seq single-end "
        psv['name'] += "SE"
    psv['title']   += psv['experiment']+" - "+psv['rep_tech'] + " (library '"+psv['library_id']+"')"
    psv['subTitle'] = psv['organism']+", "+psv['gender']+" and annotation '"+psv['annotation']+"'."
    psv['name']    += "_"+psv['experiment']+"_"+psv['rep_tech']

    # Non-file app inputs
    psv['rootR1'] = psv['experiment'] + psv['rep_tech'] + '_concatR1'
    psv['rootR2'] = psv['experiment'] + psv['rep_tech'] + '_concatR2'

    # Default locations (with adjustments)
    psv['refLoc'] = args.refLoc
    if psv['refLoc'] == REF_FOLDER_DEFAULT:
        psv['refLoc'] = REF_FOLDER_DEFAULT + '/' + psv['organism']
    psv['resultsLoc'] = args.resultsLoc
    if psv['resultsLoc'] == RESULT_FOLDER_DEFAULT:
        if psv['organism'] == 'mm10':
            psv['resultsLoc'] = RESULT_FOLDER_DEFAULT + '/' + psv['organism'] + '/' + psv['annotation']
        else:
            psv['resultsLoc'] = RESULT_FOLDER_DEFAULT + '/' + psv['organism']
    psv['resultsFolder'] = psv['resultsLoc'] + '/' + psv['experiment'] + '/' + psv['rep_tech']

    return psv




#######################
def main():

    args = get_args()

    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' %(args.experiment)
    response = dxencode.encoded_get(url, AUTHID, AUTHPW)
    exp = response.json()

    if not exp.get('replicates') or len(exp['replicates']) < 1:
        print "No replicates found in %s\n%s" % ( args.experiment, exp )
        sys.exit(1)

    replicate = "rep%s_%s" % (args.br, args.tr)

    reps_mapping = dxencode.choose_mapping_for_experiment(exp)
    # could try to do all replicates here
    try:
        mapping = reps_mapping[(args.br,args.tr)]
    except KeyError:
        print "Specified replicate: %s could not be found in mapping." % replicate
        print reps_mapping
        sys.exit(1)

    mapping['replicate'] = replicate

    try:
        mapping['genome'] = GENOME_MAPPING[mapping.get('organism', "Not Found")]

    except KeyError:
        print "Organism %s not currently supported" % mapping['organism']
        sys.exit(1)

    if mapping['unpaired'] and not mapping['paired']:
        pairedEnd = False
    elif mapping['paired'] and not mapping['unpaired']:
        pairedEnd = True
    elif not mapping['unpaired'] and not mapping['paired']:
        print "Replicate has no reads either paired or unpaired"
        print mapping
        sys.exit(1)
    else:
        print "Replicate has both paired(%s) and unpaired(%s) reads, quitting." % (len(mapping['paired'], len(mapping['unpaired'])))
        print mapping
        sys.exit(1)

    psv = pipeline_specific_vars(args, mapping, pairedEnd)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()


    ## TODO this is a bunch of ugly
    if pairedEnd:
        paired_fqs = {
            '1': [],
            '2': []
        }
        read1s = []
        read2s = []
        for (p1, p2) in mapping['paired']:
            paired_fqs[p1['paired_end']].append(p1['accession']+".fastq.gz")
            paired_fqs[p2['paired_end']].append(p2['accession']+".fastq.gz")
            read1s.append(p1['accession'])
            read2s.append(p2['accession'])
        pipePath = STEP_ORDER['pe']
        print "Generating workflow steps (paired-end)..."
    else:
        unpaired_fqs = [ f['accession']+".fastq.gz" for f in mapping['unpaired'] ]
        pipePath = STEP_ORDER['se']

    pipeSteps = STEPS
    file_globs = FILE_GLOBS

    print "Checking for prior results..."

    priors = dxencode.find_prior_results(pipePath,pipeSteps,psv['resultsFolder'],file_globs, projectId)

    if pairedEnd:
        priors['reads1'] = dxencode.find_file_set(paired_fqs["1"], projectId)
        priors['reads2'] = dxencode.find_file_set(paired_fqs["2"], projectId)
        submitted = {
            'reads1': read1s,
            'reads2': read2s
        }
    else:
        priors['reads1'] = dxencode.find_file_set(unpaired_fqs, projectId)
        submitted = {
            'reads1': mapping['unpaired'],
        }


    print "Determining steps to run..."
    #print priors
    #sys.exit(1)
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(pipePath,pipeSteps, priors, deprecateFiles, projectId, verbose=True)

    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId, verbose=True)

    if len(stepsToDo):
        print "Pipeline incomplete, please resubmit jobs: %s" % stepsToDo
        sys.exit(0)

    to_submit = [ k for k in priors.keys() if k.find('read') < 0 and POST_TEMPLATES.get(k) ]
    n = 0 # skip reads
    print "Attempting to submit %s files to args.experiment" % len(to_submit)
    while(to_submit):
        if n > len(priors) * len(priors):
            print "Too many itereations: %s" % priors
            break
        token = to_submit.pop(0)
        print "%s %s" % (token, priors[token])
        f_ob = POST_TEMPLATES.get(token, None)
        if f_ob:
            derive_check = f_ob.get('derived_from', [])
            if derive_check:
                derived = [ submitted[f] for f in derive_check if submitted.get(f) ]
                if not derived:
                    to_submit.append(token)
                    continue
                else:
                    f_ob['derived_from'] = list(itertools.chain(*derived))
            dxFile = dxpy.DXFile(dxid=priors[token])
            print "Post File: %s %s" % (token, dxFile.name)
            f_ob['dataset'] = args.experiment
            f_ob['lab'] = exp['lab']['@id']
            f_ob['award'] = exp['award']['@id']
            f_ob['assembly'] = mapping['genome']
            f_ob['annotation'] = args.annotation
            f_ob['notes'] = json.dumps(dxencode.create_notes(dxFile, dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)')))
            print json.dumps(f_ob, sort_keys=True, indent=4, separators=(',',': '))
            if args.test:
                fake_acc = 'ENCFF%03dAAA' % n
                print "Fake submission: %s" % fake_acc
                submitted[token] = [ fake_acc ]
            n += 1

    # Exit if test only
    if args.test:
        print "Fake submitted %s files." % n
    if args.test:
        sys.exit(0)


if __name__ == '__main__':
    main()

