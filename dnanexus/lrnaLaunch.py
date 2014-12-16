#!/usr/bin/env python
import argparse
import os
import sys
#import subprocess
#from datetime import datetime
#import json

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

    ap.add_argument('-r', '--replicate',
                    help="Replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=True)

    ap.add_argument('-tr', '--techrep',
                    help="Technical replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=False)

    ap.add_argument('-1', '--reads1',
                    help='Only reads fastq file or first of pair-end reads.',
                    nargs='+',
                    required=True)

    ap.add_argument('-2', '--reads2',
                    help='The second of paired-end reads files.',
                    nargs='+',
                    required=False)

    ### PIPELINE SPECIFIC
    ap.add_argument('-l', '--library',
                    help='ENCODE accession of biosample library (for BAM header)',
                    required=True)
    ### PIPELINE SPECIFIC

    ap.add_argument('-o', '--organism',
                    help="Organism to map to (default: '" + GENOME_DEFAULT + "')",
                    #choices=['hg19', 'hg38', 'mm10'],
                    choices=['hg19','mm10'],
                    default=GENOME_DEFAULT,
                    required=False)

    ap.add_argument('-g', '--gender',
                    help="Gender of sample (default: 'male')",
                    choices=['male', 'female'],
                    default='male',
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

    ap.add_argument('--run',
                    help='Run the workflow after assembling it.',
                    action='store_true',
                    required=False)

    ap.add_argument('--test',
                    help='Test run only, do not launch anything.',
                    action='store_true',
                    required=False)

    ap.add_argument('--force',
                    help='Force rerunning all steps.',
                    action='store_true',
                    required=False)

    #ap.add_argument('-x', '--export',
    #                help='Export generic Workflow (no inputs) to DNA Nexus project',
    #                action='store_true',
    #                required=False)

    return ap.parse_args()

def pipeline_specific_vars(args, pairedEnd):
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
    psv['organism']   = args.organism
    psv['gender']     = args.gender
    psv['annotation'] = args.annotation
    if psv['organism'] == 'hg19' and psv['annotation'] not in [ANNO_DEFAULT]:
        print psv['organism']+" has no "+psv['annotation']+" annotation."
        sys.exit(1)
    if psv['organism'] == 'mm10' and psv['annotation'] not in ['M2','M3','M4']:
        print psv['organism']+" has no '"+psv['annotation']+"' annotation."
        sys.exit(1)
    psv['experiment'] = args.experiment
    psv['replicate']  = str(args.replicate)
    psv['rep_tech']   = 'rep' + str(args.replicate) + '_' + str(args.techrep)
    psv['library_id'] = args.library
    psv['nthreads']   = 8
    psv['rnd_seed']   = 12345
    psv['paired']  = pairedEnd

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


def find_ref_files(priors,psv):
    '''Locates all reference files based upon gender, organism and annotation.'''
    refFiles = {}
    topIx = psv['refLoc']+'/'+GENOME_REFERENCES['tophat_index'][psv['organism']][psv['gender']][psv['annotation']]
    topIxFid = dxencode.find_file(topIx,REF_PROJECT_DEFAULT)
    if topIxFid == None:
        sys.exit("ERROR: Unable to locate TopHat index file '" + topIx + "'")
    else:
        priors['tophat_index'] = topIxFid

    starIx = psv['refLoc']+'/'+GENOME_REFERENCES['star_index'][psv['organism']][psv['gender']][psv['annotation']]
    starIxFid = dxencode.find_file(starIx,REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    rsemIx = psv['refLoc']+'/'+GENOME_REFERENCES['rsem_index'][psv['organism']][psv['annotation']]
    rsemIxFid = dxencode.find_file(rsemIx,REF_PROJECT_DEFAULT)
    if rsemIxFid == None:
        sys.exit("ERROR: Unable to locate RSEM index file '" + rsemIx + "'")
    else:
        priors['rsem_index'] = rsemIxFid

    chromSizes = psv['refLoc']+'/'+GENOME_REFERENCES['chrom_sizes'][psv['organism']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chrom_sizes'] = chromSizesFid


#######################
def main():

    args = get_args()
    if len(args.reads1) < 1:
        sys.exit('Need to have at at least 1 reads1 fastq file.')
    if args.reads2 == None:
        args.reads2 = []  # Normalize
    pairedEnd = False
    if len(args.reads2) != 0:
        pairedEnd = True
    psv = pipeline_specific_vars(args, pairedEnd)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()

    #print "Building apps dictionary..."
    pipePath = STEP_ORDER['se']
    if pairedEnd:
        pipePath = STEP_ORDER['pe']
    #pipeSteps, file_globs = dxencode.build_simple_steps(pipePath,projectId,verbose=True)
    pipeSteps = STEPS
    file_globs = FILE_GLOBS

    print "Checking for prior results..."
    # Check if there are previous results
    # Perhaps reads files are already there?
    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    if not args.test:
        if not dxencode.project_has_folder(project, psv['resultsFolder']):
            project.new_folder(psv['resultsFolder'],parents=True)

    priors = dxencode.find_prior_results(pipePath,pipeSteps,psv['resultsFolder'],file_globs, projectId)

    print "Checking for read files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    inputs = {}
    inputs['Reads1'] = dxencode.find_and_copy_read_files(priors, args.reads1, args.test, 'reads1', \
                                                            psv['resultsFolder'], False, projectId)
    inputs['Reads2'] = dxencode.find_and_copy_read_files(priors, args.reads2, args.test, 'reads2', \
                                                            psv['resultsFolder'], False, projectId)

    print "Looking for reference files..."
    find_ref_files(priors,psv)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(pipePath,pipeSteps, priors, deprecateFiles, projectId, \
                                                                                force=args.force)

    # Report the plans
    dxencode.report_plans(psv, inputs, GENOME_REFERENCES.keys(), deprecateFiles, priors, \
                                                                pipePath, stepsToDo, pipeSteps)
    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId, verbose=True)

    if len(deprecateFiles) > 0 and not args.test:
        oldFolder = psv['resultsFolder']+"/deprecated"
        print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '"+oldFolder+"/'..."
        dxencode.move_files(deprecateFiles,oldFolder,projectId)

    if args.test:
        print "Testing workflow assembly..."
    else:
        print "Assembling workflow..."
    wf = dxencode.create_workflow(stepsToDo, pipeSteps, priors, psv, projectId, test=args.test)

    # Exit if test only
    if args.test:
        print "TEST ONLY - exiting."
        sys.exit(0)

    # Roll out to pad and possibly launch
    dxencode.launchPad(wf,projectId,psv,args.run)

    print "(success)"

if __name__ == '__main__':
    main()

