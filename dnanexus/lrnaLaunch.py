#!/usr/bin/env python
import argparse
import os
import sys
#import subprocess
#from datetime import datetime
#import json

import dxpy
import dxencode as dxencode
#from dxencode import dxencode as dxencode

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
    'se': [ 'concatR1',             'tophatSe', 'topBwSe', 'starSe', 'starBwSe', 'rsem' ],
    'pe': [ 'concatR1', 'concatR2', 'tophatPe', 'topBwPe', 'starPe', 'starBwPe', 'rsem' ]
    }

STEPS = {
    # for each step: app, list of any params, inputs and results (both as fileToken: app_obj_name)
    # TODO: Any results files not listed here would not be 'deprecated' on reruns.
    'concatR1': {
                'app':     'concat-fastqs',
                'params':  { 'rootR1': 'outfile_root' },
                'inputs':  { 'reads1_set': 'reads_set' },
                'results': { 'reads1': 'reads' }
                },
    'concatR2': {
                'app':     'concat-fastqs',
                'params':  { 'rootR2': 'outfile_root' },
                'inputs':  { 'reads2_set': 'reads_set' },
                'results': { 'reads2': 'reads' }
                },
    'tophatSe': {
                'app':     'align-tophat-se',
                'params':  { 'library': 'library_id' }, #, 'nthreads'
                'inputs':  { 'reads1': 'reads', 'tophatIndex': 'tophat_index' },
                'results': { 'topBam': 'genome_bam' }
                },
    'tophatPe': {
                'app':     'align-tophat-pe',
                'params':  { 'library': 'library_id' }, #, 'nthreads'
                'inputs':  { 'reads1': 'reads_1', 'reads2': 'reads_2', 'tophatIndex': 'tophat_index' },
                'results': { 'topBam': 'genome_bam' }
                },
    'topBwSe':  {
                'app':     'bam-to-bigwig-unstranded',
                'inputs':  { 'topBam': 'bam_file', 'chromSizes': 'chrom_sizes' },
                'results': { 'topBwAll': 'all_unstranded_bw', 'topBwUniq': 'unique_unstranded_bw' }
                },
    'topBwPe':  {
                'app':     'bam-to-bigwig-stranded',
                'inputs':  { 'topBam': 'bam_file', 'chromSizes': 'chrom_sizes' },
                'results': { 'topBwMinusAll': 'all_minus_bw', 'topBwMinusUniq': 'uniq_minus_bw',
                              'topBwPlusAll': 'all_plus_bw',   'topBwPlusUniq': 'all_plus_bw' }
                },
    'starSe':   {
                'app':     'align-star-se',
                'params':  { 'library': 'library_id' }, #, 'nthreads'
                'inputs':  { 'reads1': 'reads', 'starIndex': 'star_index' },
                'results': { 'starGenoBam': 'genome_bam', 'starAnnoBam': 'annotation_bam',
                                                                'starLog': 'star_log' }
                },
    'starPe':   {
                'app':     'align-star-pe',
                'params':  { 'library': 'library_id' }, #, 'nthreads'
                'inputs':  { 'reads1': 'reads_1', 'reads2': 'reads_2', 'starIndex': 'star_index' },
                'results': { 'starGenoBam': 'genome_bam', 'starAnnoBam': 'annotation_bam',
                                                                'starLog': 'star_log' }
                },
    'starBwSe': {
                'app':     'bam-to-bigwig-unstranded',
                'inputs':  { 'starGenoBam': 'bam_file', 'chromSizes': 'chrom_sizes' },
                'results': { 'starBwAll': 'all_unstranded_bw', 'starBwUniq': 'unique_unstranded_bw' }
                },
    'starBwPe': {
                'app':     'bam-to-bigwig-stranded',
                'inputs':  { 'starGenoBam': 'bam_file', 'chromSizes': 'chrom_sizes' },
                'results': { 'starBwMinusAll': 'all_minus_bw','starBwMinusUniq': 'uniq_minus_bw',
                              'starBwPlusAll': 'all_plus_bw',  'starBwPlusUniq': 'uniq_plus_bw' }
                },
    'rsem':     {
                'app':     'quant-rsem',
                'params':  { 'pairedEnd': 'paired' },  #, 'nthreads', 'rnd_seed'
                'inputs':  { 'starAnnoBam': 'annotation_bam', 'rsemIndex': 'rsem_index' },
                'results': { 'rsemIso': 'transcript_quant', 'rsemGene': 'genomic_quant' }
                }
    }

FILE_GLOBS = {
    # For looking up previous result files, use wild-cards
    'reads1':          "/*_concatR1.fq.gz",
    'reads2':          "/*_concatR2.fq.gz",
    'topBam':          "/*_tophat.bam",
    'topBwMinusAll':   "/*_tophat_minusAll.bw",
    'topBwMinusUniq':  "/*_tophat_minusUniq.bw",
    'topBwPlusAll':    "/*_tophat_plusAll.bw",
    'topBwPlusUniq':   "/*_tophat_plusUniq.bw",
    'topBwAll':        "/*_tophat_all.bw",
    'topBwUniq':       "/*_tophat_uniq.bw",
    'starGenoBam':     "/*_star_genome.bam",
    'starAnnoBam':     "/*_star_anno*.bam",
    'starLog':         "/*_Log.final.out",
    'starBwMinusAll':  "/*_star_genome_minusAll.bw",
    'starBwMinusUniq': "/*_star_genome_minusUniq.bw",
    'starBwPlusAll':   "/*_star_genome_plusAll.bw",
    'starBwPlusUniq':  "/*_star_genome_plusUniq.bw",
    'starBwAll':       "/*_star_genome_all.bw",
    'starBwUniq':      "/*_star_genome_uniq.bw",
    'rsemIso':         "/*_rsem.isoforms.results",
    'rsemGene':        "/*_rsem.genes.results"
    }

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should remove annotation if only one per genome
    # TODO: should use ACCESSION based fileNames
    'tophatIndex':  {
                    'hg19': {
                            'female':   {
                                        'v19': 'hg19_female_v19_ERCC_tophatIndex.tgz'
                                        },
                            'male':     {
                                        'v19': 'hg19_male_v19_ERCC_tophatIndex.tgz'
                                        }
                            },
                    'mm10': {
                            'female':   {
                                        'M2':  'mm10_female_M2_ERCC_tophatIndex.tgz',
                                        'M3':  'mm10_female_M3_ERCC_tophatIndex.tgz',
                                        'M4':  'mm10_female_M4_ERCC_tophatIndex.tgz'
                                        },
                            'male':     {
                                        'M2':  'mm10_male_M2_ERCC_tophatIndex.tgz',
                                        'M3':  'mm10_male_M3_ERCC_tophatIndex.tgz',
                                        'M4':  'mm10_male_M4_ERCC_tophatIndex.tgz'
                                        }
                            }
                    },
    'starIndex':    {
                    'hg19': {
                            'female':   {
                                        'v19': 'hg19_female_v19_ERCC_starIndex.tgz'
                                        },
                            'male':     {
                                        'v19': 'hg19_male_v19_ERCC_starIndex.tgz'
                                        }
                            },
                    'mm10': {
                            'female':   {
                                        'M2':  'mm10_female_M2_ERCC_starIndex.tgz',
                                        'M3':  'mm10_female_M3_ERCC_starIndex.tgz',
                                        'M4':  'mm10_female_M4_ERCC_starIndex.tgz'
                                        },
                            'male':     {
                                        'M2':  'mm10_male_M2_ERCC_starIndex.tgz',
                                        'M3':  'mm10_male_M3_ERCC_starIndex.tgz',
                                        'M4':  'mm10_male_M4_ERCC_starIndex.tgz'
                                        }
                            }
                    },
    'rsemIndex':    {
                    'hg19': {
                            'v19': 'hg19_male_v19_ERCC_rsemIndex.tgz'
                            },
                    'mm10': {
                            'M2':  'mm10_male_M2_ERCC_rsemIndex.tgz',
                            'M3':  'mm10_male_M3_ERCC_rsemIndex.tgz',
                            'M4':  'mm10_male_M4_ERCC_rsemIndex.tgz'
                            }
                    },
    'chromSizes':   {
                    'hg19': {
                            'female':   'female.hg19.chrom.sizes',
                            'male':     'male.hg19.chrom.sizes'
                            },
                    'mm10': {
                            'female':   'female.mm10.chrom.sizes',
                            'male':     'male.mm10.chrom.sizes'
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
    psv['library']    = args.library
    psv['pairedEnd']  = pairedEnd

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
    psv['title']   += psv['experiment']+" - "+psv['rep_tech'] + " (library '"+psv['library']+"')"
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
    topIx = psv['refLoc']+'/'+GENOME_REFERENCES['tophatIndex'][psv['organism']][psv['gender']][psv['annotation']]
    topIxFid = dxencode.find_file(topIx,REF_PROJECT_DEFAULT)
    if topIxFid == None:
        sys.exit("ERROR: Unable to locate TopHat index file '" + topIx + "'")
    else:
        priors['tophatIndex'] = topIxFid

    starIx = psv['refLoc']+'/'+GENOME_REFERENCES['starIndex'][psv['organism']][psv['gender']][psv['annotation']]
    starIxFid = dxencode.find_file(starIx,REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['starIndex'] = starIxFid

    rsemIx = psv['refLoc']+'/'+GENOME_REFERENCES['rsemIndex'][psv['organism']][psv['annotation']]
    rsemIxFid = dxencode.find_file(rsemIx,REF_PROJECT_DEFAULT)
    if rsemIxFid == None:
        sys.exit("ERROR: Unable to locate RSEM index file '" + rsemIx + "'")
    else:
        priors['rsemIndex'] = rsemIxFid

    chromSizes = psv['refLoc']+'/'+GENOME_REFERENCES['chromSizes'][psv['organism']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chromSizes'] = chromSizesFid


#######################
def main():

    args = get_args()
    if len(args.reads1) < 1:
        sys.exit('Need to have at least 1 replicate file.')
    if args.reads2 == None:
        args.reads2 = []  # Normalize
    pairedEnd = False
    if len(args.reads2) != 0:
        pairedEnd = True
    psv = pipeline_specific_vars(args, pairedEnd)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()

    print "Checking for prior results..."
    # Check if there are previous results
    # Perhaps reads files are already there?
    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    pipePath = STEP_ORDER['se']
    if pairedEnd:
        pipePath = STEP_ORDER['pe']
    if not args.test:
        if not dxencode.project_has_folder(project, psv['resultsFolder']):
            project.new_folder(psv['resultsFolder'],parents=True)
    priors = dxencode.find_prior_results(pipePath,STEPS,psv['resultsFolder'],FILE_GLOBS, projectId)
    
    print "Checking for read files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    reads1 = dxencode.find_and_copy_read_files(priors, args.reads1, args.test, 'reads1', \
                                                            psv['resultsFolder'], False, projectId)
    reads2 = dxencode.find_and_copy_read_files(priors, args.reads2, args.test, 'reads2', \
                                                            psv['resultsFolder'], False, projectId)

    print "Looking for reference files..."
    find_ref_files(priors,psv)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(pipePath,STEPS, priors, deprecateFiles, projectId, \
                                                                                force=args.force)

    # Report the plans
    dxencode.report_plans(psv, reads1, reads2, GENOME_REFERENCES.keys(), deprecateFiles, priors, \
                                                                pipePath, stepsToDo, STEPS)
    print "Running '"+psv['title']+"'"
    print "     on "+psv['subTitle']
    if pairedEnd:
        print "- Reads1: "
    else:
        print "- Reads: "
    for fid in reads1:
        print "  " + dxencode.file_path_from_fid(fid)
    if pairedEnd:
        print "- Reads2: "
        for fid in reads2:
            print "  " + dxencode.file_path_from_fid(fid)
    print "- Reference files:"
    for token in GENOME_REFERENCES.keys():
        print "  " + dxencode.file_path_from_fid(priors[token],True)
    print "- Results written to: " + args.project + ":" +psv['resultsFolder'] +'/'
    if len(stepsToDo) == 0:
        print "* All expected results are in the results folder, so there is nothing to do."
        print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
        print "  rerun all steps; or remove suspect results from the folder before launching."
        sys.exit(0)
    else:
        print "- Steps to run:"
        for step in pipePath:
            if step in stepsToDo:
                print "  * "+STEPS[step]['app']+" will be run"
            else:
                if not step.find('concat') == 0:
                    print "    "+STEPS[step]['app']+" has already been run"

    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId, verbose=True)

    if len(deprecateFiles) > 0:
        oldFolder = psv['resultsFolder']+"/deprecated"
        if args.test:
            print "Would move "+str(len(deprecateFiles))+" prior result file(s) to '"+oldFolder+"/'."
            for fid in deprecateFiles:
                print "  " + dxencode.file_path_from_fid(fid)
        else:
            print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '"+oldFolder+"/'..."
            dxencode.move_files(deprecateFiles,oldFolder,projectId)

    if args.test:
        print "Testing workflow assembly..."
    else:
        print "Assembling workflow..."
    wf = dxencode.create_workflow(stepsToDo, STEPS, priors, psv, projectId, test=args.test)

    # Exit if test only
    if args.test:
        print "TEST ONLY - exiting."
        sys.exit(0)
        
    # Roll out to pad and possibly launch
    dxencode.launchPad(wf,projectId,psv,args.run)        
            
    print "(success)"

if __name__ == '__main__':
    main()

