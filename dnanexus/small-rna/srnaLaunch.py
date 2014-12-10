#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
#from datetime import datetime

import dxpy
import dxencode as dxencode
#from dxencode import dxencode as dxencode
import json

# NOTES: This command-line utility will run the short RNA-seq pipeline for a single replicate
#      - All results will be written to a folder /srna/<expId>/rep<#>.
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

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''

RESULT_FOLDER_DEFAULT = '/srna'
''' This the default location to place results folders for each experiment.'''

STEP_ORDER = [ 'concat-fastqs', 'small-rna-align', 'small-rna-signals' ]

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should use ACCESSION based fileNames
    'star_index':   {
                    'hg19': {
                            'female':   'hg19_female_sRNA_starIndex.tgz',
                            'male':     'hg19_male_sRNA_starIndex.tgz'
                            },
                    'mm10': {
                            'female':   'mm10_female_sRNA_starIndex.tgz',
                            'male':     'mm10_male_sRNA_starIndex.tgz'
                            }
                    },
    'chrom_sizes':  {
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
    ap = argparse.ArgumentParser(description="Launches short RNA-seq pipeline analysis for " +
                "one replicate. Can be run repeatedly and will launch only the steps that " +
                "are needed to finish the pipeline. All results will be placed in the " +
                "folder /<resultsLoc>/<experiment>/<replicate>.")
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

    ap.add_argument('-1', '--reads',
                    help='Only reads fastq file or files.',
                    nargs='+',
                    required=True)

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

def pipeline_specific_vars(args):
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
    psv['experiment'] = args.experiment
    psv['replicate']  = str(args.replicate)
    psv['rep_tech']   = 'rep' + str(args.replicate) + '_' + str(args.techrep)
    psv['library_id'] = args.library
    psv['nthreads'] = 8
    # workflow labeling
    genderToken = "XY"
    if psv['gender'] == 'female':
        genderToken = "XX"
    psv['description'] = "The ENCODE RNA Seq pipeline for short RNA"
    psv['title'] = "short RNA-seq " + psv['experiment'] + " - "+psv['rep_tech'] + \
                   " (library '"+psv['library_id']+"') on " + psv['organism']+" - "+psv['gender']
    psv['name'] = "srna_"+psv['organism']+genderToken+"_"+psv['experiment']+"_"+psv['rep_tech']

    # Non-file app inputs
    psv['rootR1'] = psv['experiment'] + psv['rep_tech'] + '_concatReads'

    # Default locations (with adjustments)
    psv['refLoc'] = args.refLoc
    if psv['refLoc'] == REF_FOLDER_DEFAULT:
        psv['refLoc'] = REF_FOLDER_DEFAULT + '/' + psv['organism']
    psv['resultsLoc'] = args.resultsLoc
    if psv['resultsLoc'] == RESULT_FOLDER_DEFAULT:
        psv['resultsLoc'] = RESULT_FOLDER_DEFAULT +'/'+psv['organism']
    psv['resultsFolder'] = psv['resultsLoc'] + '/' + psv['experiment'] + '/' + psv['rep_tech']

    return psv


def find_ref_files(priors,psv):
    '''Locates all reference files based upon organism and gender.'''
    refFiles = {}
    starIx = psv['refLoc']+'/'+GENOME_REFERENCES['star_index'][psv['organism']][psv['gender']]
    starIxFid = dxencode.find_file(starIx,REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    chromSizes = psv['refLoc']+'/'+GENOME_REFERENCES['chrom_sizes'][psv['organism']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chrom_sizes'] = chromSizesFid

#######################
def main():

    args = get_args()
    if len(args.reads) < 1:
        sys.exit('Need to have at least 1 replicate file.')
    psv = pipeline_specific_vars(args)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()

    print "Building apps dictionary..."
    #pipeSteps = STEPS
    #file_globs = FILE_GLOBS
    pipeSteps, file_globs = dxencode.build_simple_steps(STEP_ORDER,projectId)
    

    print "Checking for prior results..."
    # Check if there are previous results
    # Perhaps reads files are already there?
    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    if not args.test:
        if not dxencode.project_has_folder(project, psv['resultsFolder']):
            project.new_folder(psv['resultsFolder'],parents=True)
    priors = dxencode.find_prior_results(STEP_ORDER,pipeSteps,psv['resultsFolder'],file_globs,projectId)

    print "Checking for read files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    reads = dxencode.find_and_copy_read_files(priors, args.reads, args.test, 'reads', \
                                                            psv['resultsFolder'], False, projectId)

    print "Looking for reference files..."
    find_ref_files(priors,psv)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(STEP_ORDER, pipeSteps, priors, deprecateFiles, \
                                                                    projectId, force=args.force)
    # Report the plans
    dxencode.report_plans(psv, reads, None, GENOME_REFERENCES.keys(), deprecateFiles, priors, \
                                                                STEP_ORDER, stepsToDo, pipeSteps)
    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId,verbose=True)

    if len(deprecateFiles) > 0 and not args.test:
        oldFolder = psv['resultsFolder']+"/deprecated"
        print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '"+oldFolder+"/'..."
        exencode.move_files(deprecateFiles,oldFolder,projectId)

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

