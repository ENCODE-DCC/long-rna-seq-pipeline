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

GENOMES_SUPPORTED = ['hg19', 'mm10']
GENOME_DEFAULT = 'hg19'
''' This the default Genome that long RNA-seq experiments are mapped to.'''

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

RESULT_FOLDER_DEFAULT = '/srna/'
''' This the default location to place results folders for each experiment.'''

REP_STEP_ORDER = [ "concat-fastqs", "small-rna-align", "small-rna-signals" ]
'''The (artifically) linear order of all pipeline steps.'''

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should use ACCESSION based fileNames
    "star_index":   {
                    "hg19": {
                            "female":   "hg19_female_sRNA_starIndex.tgz",
                            "male":     "hg19_male_sRNA_starIndex.tgz"
                            },
                    "mm10": {
                            "female":   "mm10_female_sRNA_starIndex.tgz",
                            "male":     "mm10_male_sRNA_starIndex.tgz"
                            }
                    },
    "chrom_sizes":  {
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
    ap = argparse.ArgumentParser(description="Launches short RNA-seq pipeline analysis for " +
                "one replicate. Can be run repeatedly and will launch only the steps that " +
                "are needed to finish the pipeline. All results will be placed in the " +
                "folder /<resultsLoc>/<experiment>/<replicate>.")
    ### PIPELINE SPECIFIC

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    required=True)

    ap.add_argument('--br', '--biological-replicate',
                    help="Biological replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=True)

    ap.add_argument('--tr', '--technical-replicate',
                    help="Technical replicate number (default: 1)",
                    type=int,
                    default='1',
                    required=False)

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    default=PROJECT_DEFAULT,
                    required=False)

    ap.add_argument('--refLoc',
                    help="The location to find reference files (default: '" + \
                            dxencode.REF_PROJECT_DEFAULT + ":" + dxencode.REF_FOLDER_DEFAULT + "')",
                    default=dxencode.REF_FOLDER_DEFAULT,
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

def pipeline_specific_vars(args,verbose=False):
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

    # Start with dict containing common variables
    print "Retrieving experiment specifics..."
    psv = dxencode.common_variables(args,RESULT_FOLDER_DEFAULT,controls=False)
    
    # Now add pipline specific variables and tests
    
    # Paired ends?
    if psv['paired_end']:
        print "Small-RNA is always expected to be single-end but mapping says otherwise."
        print mapping
        sys.exit(1)

    # Non-file app inputs
    # Ugly fixup of standard names in other pipelines.
    for rep in psv['reps'].values():
        del rep['rootR1']
        rep['outfile_root'] = psv['experiment'] + rep['rep_tech'] + '_concatReads'

    # Some specific settings
    psv['nthreads']   = 8

    # run will either be for combined or single rep.
    if not psv['combined']:
        run = psv['reps']['a']  # If not combined then run will be for the first (only) replicate
    else:
        run = psv
        print "Small-RNA-seq pipeline currently does not support combined-replicate processing."
        print mapping
        sys.exit(1)
        
    # workflow labeling
    psv['description'] = "The ENCODE RNA Seq pipeline for short RNA"
    genderToken = "XY"
    if psv['gender'] == 'female':
        genderToken = "XX"
    run['title'] = "short RNA-seq " + psv['experiment'] + " - "+run['rep_tech'] + \
                   " (library '"+run['library_id']+"') on " + psv['genome']+" - "+psv['gender']
    run['name'] = "srna_"+psv['genome']+genderToken+"_"+psv['experiment']+"_"+run['rep_tech']

    if verbose:
        print "Pipeline Specific Vars:"
        print json.dumps(psv,indent=4)
    return psv


def find_ref_files(priors,psv):
    '''Locates all reference files based upon organism and gender.'''
    refFiles = {}
    starIx = psv['refLoc']+GENOME_REFERENCES['star_index'][psv['genome']][psv['gender']]
    starIxFid = dxencode.find_file(starIx,dxencode.REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    chromSizes = psv['refLoc']+GENOME_REFERENCES['chrom_sizes'][psv['genome']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,dxencode.REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chrom_sizes'] = chromSizesFid
    psv['ref_files'] = GENOME_REFERENCES.keys()

#######################
def main():

    args = get_args()
    print "Retrieving pipeline specifics..."
    psv = pipeline_specific_vars(args)
    
    project = dxencode.get_project(args.project)
    projectId = project.get_id()

    print "Building apps dictionary..."
    pipeRepSteps, file_globs = dxencode.build_simple_steps(REP_STEP_ORDER,projectId)
    for rep in psv['reps'].values():
        rep['path'] = REP_STEP_ORDER

    # finding fastqs and prior results in a stadardized way
    dxencode.finding_rep_inputs_and_priors(psv,pipeRepSteps,file_globs,project,args.test)

    # finding pipeline specific reference files in a stadardized way
    dxencode.find_all_ref_files(psv,find_ref_files)

    # deterine steps to run in a stadardized way
    dxencode.determine_steps_needed(psv,pipeRepSteps,None,projectId, args.force)

    # Preperation is done. At this point on we either run rep 'a' or combined.
    run = psv['reps']['a']
    run['steps'] = pipeRepSteps
        
    dxencode.report_build_launch(psv, run, projectId, test=args.test, launch=args.run)
            
    print "(success)"

if __name__ == '__main__':
    main()

