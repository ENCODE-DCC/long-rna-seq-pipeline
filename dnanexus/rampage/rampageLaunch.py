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
GENOME_DEFAULTS = { 'human': 'hg19', 'mouse': 'mm10' }
GENOME_DEFAULT = 'hg19'
''' This the default Genome that long RNA-seq experiments are mapped to.'''

ANNO_DEFAULTS = {'hg19': 'v19', 'mm10':'M4' }
ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                 'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
ANNO_DEFAULT = ANNO_DEFAULTS['hg19']
''' Multiple annotations might be supported for each genome.'''

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''

RESULT_FOLDER_DEFAULT = '/rampage/'
''' This the default location to place results folders for each experiment.'''

# NOTE '/lrna/' or '/run/' is safer and more efficient than '/', but also more brittle
#CONTROL_ROOT_FOLDER = '/'
CONTROL_ROOT_FOLDER = '/lrna/'
''' Rampage requires a control file which may be discoverable.'''
CONTROL_FILE_GLOB = '*_star_genome.bam'


STEP_ORDER = [ "concatR1", "concatR2", "rampage-align-pe", "rampage-signals", "rampage-peaks" ]
#STEP_ORDER = [ "rampage-align-pe", "rampage-signals", "rampage-peaks" ]
'''The (artifically) linear order of all pipeline steps.'''

STEPS = {
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
    "rampage-align-pe": {
        "app":     "rampage-align-pe",
        "params":  { "library_id": "library_id", "nthreads": "nthreads" },
        "inputs":  { "reads1": "reads_1", 
                     "reads2": "reads_2", 
                     "star_index": "star_index" },
        "results": { "rampage_marked_bam": "rampage_marked_bam", 
                     "rampage_star_log": "rampage_star_log" }
    },
    "rampage-signals": {
        "app":    "rampage-signals",
        "params":  {},
        "inputs":  { "chrom_sizes": "chrom_sizes", 
                     "rampage_marked_bam": "rampage_marked_bam" },
        "results": { "all_plus_bw":    "all_plus_bw",    
                     "all_minus_bw":    "all_minus_bw",
                     "unique_plus_bw": "unique_plus_bw", 
                     "unique_minus_bw": "unique_minus_bw",
                     "all_plus_bg":    "all_plus_bg",    
                     "all_minus_bg":    "all_minus_bg", 
                     "unique_plus_bg": "unique_plus_bg", 
                     "unique_minus_bg": "unique_minus_bg" }
    },
    "rampage-peaks": {
        "app":     "rampage-peaks",
        "params":  { "nthreads": "nthreads" },
        "inputs":  { "control_bam": "control_bam", 
                     "gene_annotation": "gene_annotation", 
                     "chrom_sizes": "chrom_sizes",
                     "rampage_marked_bam": "rampage_marked_bam" },
        "results": { "rampage_peaks_bed": "rampage_peaks_bed",
                     "rampage_peaks_bb": "rampage_peaks_bb",
                     "rampage_peaks_gtf": "rampage_peaks_gtf" }
    }
}

FILE_GLOBS = {
    "reads1":             "/*_concatR1.fq.gz",
    "reads2":             "/*_concatR2.fq.gz",
    "all_plus_bw":        "/*_rampage_5p_plusAll.bw",
    "rampage_marked_bam": "/*_rampage_star_marked.bam",
    "all_minus_bg":       "/*_rampage_5p_minusAll.bg",
    "unique_plus_bg":     "/*_rampage_5p_plusUniq.bg",
    "rampage_peaks_bed":  "/*_rampage_peaks.bed",
    "rampage_peaks_gtf":  "/*_rampage_peaks.gtf",
    "all_plus_bg":        "/*_rampage_5p_plusAll.bg",
    "rampage_star_log":   "/*_rampage_star_Log.final.out",
    "all_minus_bw":       "/*_rampage_5p_minusAll.bw",
    "rampage_peaks_bb":   "/*_rampage_peaks.bb",
    "unique_plus_bw":     "/*_rampage_5p_plusUniq.bw",
    "unique_minus_bg":    "/*_rampage_5p_minusUniq.bg",
    "unique_minus_bw":    "/*_rampage_5p_minusUniq.bw"
}
GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should remove annotation if only one per genome
    # TODO: should use ACCESSION based fileNames
    "star_index":   {
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
    "gene_annotation":   {
                    "hg19": { "v19": "gencode.v19.annotation.gtf.gz" },
                    "mm10": {
                              "M2":  "gencode.vM2.annotation.gtf.gz",
                              "M3":  "gencode.vM3.annotation.gtf.gz",
                              "M4":  "gencode.vM4.annotation.gtf.gz"
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

    ap.add_argument('-c', '--control',
                    help='The control bam for peak calling.',
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

    psv = {}
    psv['experiment'] = args.experiment
    psv['biorep']  = str(args.br)
    psv['rep_tech']   = 'rep' + str(args.br) + '_' + str(args.tr)

    mapping = dxencode.get_mapping(args.experiment,args.br,args.tr)

    # Only supported genomes
    if mapping['organism'] in GENOME_DEFAULTS:
        psv['genome'] = GENOME_DEFAULTS[mapping['organism']]
    else:
        print "Organism %s not currently supported" % mapping['organism']
        sys.exit(1)

    # Could be multiple annotations supported per genome
    psv['annotation'] = args.annotation
    if psv['genome'] != GENOME_DEFAULT and psv['annotation'] == ANNO_DEFAULT:
        psv['annotation'] = ANNO_DEFAULTS[psv['genome']]
    if psv['annotation'] not in ANNO_ALLOWED[psv['genome']]:
        print psv['genome']+" has no "+psv['annotation']+" annotation."
        sys.exit(1)
    
    # Paired ends?  Read files?
    psv['paired_end'] = dxencode.load_fastqs_from_mapping(psv, mapping, controls=True)
    if not psv['paired_end']:
        print "Rampage is always expected to be paired-end but mapping says otherwise."
        print mapping
        sys.exit(1)

    # And the rest
    psv['gender']= mapping['sex']
    psv['library_id'] = mapping['library']
    psv['project']    = args.project
    psv['nthreads']   = 8

    # Non-file app inputs
    psv['rootR1'] = psv['experiment'] + psv['rep_tech'] + '_concatR1'
    psv['rootR2'] = psv['experiment'] + psv['rep_tech'] + '_concatR2'

    # workflow labeling
    psv['description'] = "The ENCODE Rampage RNA pipeline for long RNAs"
    psv['name'] = "rampage_"+psv['genome']
    if psv['genome'] == 'mm10':
        psv['name'] += psv['annotation']
    if psv['gender'] == 'female':
        psv['name'] += "XX"
    else:
        psv['name'] += "XY"
    psv['title'] = "Rampage RNA "
    psv['title'] += psv['experiment'] + " - "+psv['rep_tech'] + \
                    " (library '"+psv['library_id']+"') on " + psv['genome']+" - "+psv['gender']
    psv['name']    += "_"+psv['experiment']+"_"+psv['rep_tech']

    # Default locations (with adjustments)
    psv['refLoc'] = args.refLoc
    if psv['refLoc'] == REF_FOLDER_DEFAULT:
        psv['refLoc'] = REF_FOLDER_DEFAULT + psv['genome'] + '/'
    if not psv['refLoc'].endswith('/'):
        psv['refLoc'] += '/' 
    psv['resultsLoc'] = args.resultsLoc
    if psv['resultsLoc'] == RESULT_FOLDER_DEFAULT:
        if psv['genome'] == 'mm10':
            psv['resultsLoc'] = RESULT_FOLDER_DEFAULT + psv['genome'] + '/' + psv['annotation'] + '/'
        else:
            psv['resultsLoc'] = RESULT_FOLDER_DEFAULT + psv['genome'] + '/'
    if not psv['resultsLoc'].endswith('/'):
        psv['resultsLoc'] += '/' 
    psv['resultsFolder'] = psv['resultsLoc'] + psv['experiment'] + '/' + psv['rep_tech'] + '/'

    if verbose:
        print "Pipeline Specific Vars:"
        print json.dumps(psv,indent=4)
    return psv


def find_ref_files(priors,psv):
    '''Locates all reference files based upon organism and gender.'''
    refFiles = {}
    starIx = psv['refLoc']+GENOME_REFERENCES['star_index'][psv['genome']][psv['gender']][psv['annotation']]
    starIxFid = dxencode.find_file(starIx,REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    anno = psv['refLoc']+GENOME_REFERENCES['gene_annotation'][psv['genome']][psv['annotation']]
    anno_fid = dxencode.find_file(anno,REF_PROJECT_DEFAULT)
    if anno_fid == None:
        sys.exit("ERROR: Unable to locate Gene Annotation file '" + anno + "'")
    else:
        priors['gene_annotation'] = anno_fid

    chromSizes = psv['refLoc']+GENOME_REFERENCES['chrom_sizes'][psv['genome']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chrom_sizes'] = chromSizesFid
    

def find_control_file(psv,project,default=None):
    '''Attempts to find an appropriate control file.'''
    # TODO Make more generic and move to dxencode.py when needed.
    
    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    for file_key in psv['controls']:
        url = '%s%s/?format=json&frame=embedded' % (SERVER,file_key)
        response = dxencode.encoded_get(url, AUTHID, AUTHPW)
        file_obj = response.json()
        rep_id = file_obj["replicate"]
        url = '%s%s/?format=json&frame=embedded' % (SERVER,rep_id)
        response = dxencode.encoded_get(url, AUTHID, AUTHPW)
        rep_obj = response.json()
        exp_id = rep_obj['experiment'].split('/')[2]
        rep_tech = "rep%s_%s" % \
                (rep_obj['biological_replicate_number'], rep_obj['technical_replicate_number'])
        # default by cheating
        control_root = CONTROL_ROOT_FOLDER
        if psv['project'] != PROJECT_DEFAULT:
            control_root = '/runs/'        # TODO: remove this cheat!
        path_n_glob = control_root + exp_id + '/' + rep_tech + '/' + CONTROL_FILE_GLOB
        target_folder = dxencode.find_folder(exp_id + '/' + rep_tech,project,control_root)
        #print "Target found [%s]" % target_folder
        if target_folder != None:
            path_n_glob = target_folder + '/' + CONTROL_FILE_GLOB
        fid = dxencode.find_file(path_n_glob,project.get_id(),multiple=False,recurse=False)
        if fid != None:
            return dxencode.file_path_from_fid(fid)
            
    if default != None:
        return default
    print "Unable to find control in search of %s" % psv['controls']
    sys.exit(1)
            
    

#######################
def main():

    args = get_args()
    psv = pipeline_specific_vars(args)

    project = dxencode.get_project(psv['project'])
    projectId = project.get_id()

    #print "Building apps dictionary..."
    #pipeSteps, file_globs = dxencode.build_simple_steps(STEP_ORDER,projectId,verbose=True)
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
    priors = dxencode.find_prior_results(STEP_ORDER,pipeSteps,psv['resultsFolder'],file_globs,projectId)

    print "Checking for input files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    inputs = {}
    inputs['Reads1'] = dxencode.find_and_copy_read_files(priors, psv['fastqs']['1'], args.test, \
                                              'reads1', psv['resultsFolder'], False, projectId)
    inputs['Reads2'] = dxencode.find_and_copy_read_files(priors, psv['fastqs']['2'], args.test, \
                                              'reads2', psv['resultsFolder'], False, projectId)
    print "Checking for control file..."
    control = find_control_file(psv,project,args.control)
    inputs['Control'] = dxencode.find_and_copy_read_files(priors, [ control ], args.test, \
                                              'control_bam', psv['resultsFolder'], False, projectId)

    print "Looking for reference files..."
    find_ref_files(priors,psv)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(STEP_ORDER, pipeSteps, priors, deprecateFiles, \
                                                                    projectId, force=args.force)
    # Report the plans
    dxencode.report_plans(psv, inputs, GENOME_REFERENCES.keys(), deprecateFiles, priors, \
                                                                STEP_ORDER, stepsToDo, pipeSteps)
    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId,verbose=True)

    if len(deprecateFiles) > 0 and not args.test:
        oldFolder = psv['resultsFolder']+"deprecated/"
        print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '"+oldFolder+"'..."
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

