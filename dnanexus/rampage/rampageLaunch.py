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
''' This the default Genomes that long RNA-seq experiments are mapped to.'''

ANNO_DEFAULTS = {'hg19': 'v19', 'mm10':'M4' }
ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                 'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
ANNO_DEFAULT = ANNO_DEFAULTS[GENOME_DEFAULT]
''' Multiple annotations might be supported for each genome.'''

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

RESULT_FOLDER_DEFAULT = '/rampage/'
''' This the default location to place results folders for each experiment.'''

# NOTE '/lrna/' or '/run/' is safer and more efficient than '/', but also more brittle
#CONTROL_ROOT_FOLDER = '/'
CONTROL_ROOT_FOLDER = '/lrna/'
''' Rampage requires a control file which may be discoverable.'''
CONTROL_FILE_GLOB = '*_star_genome.bam'


REP_STEP_ORDER = [ "concatR1", "concatR2", "rampage-align-pe", "rampage-signals", "rampage-peaks" ]
'''The (artifically) linear order of all pipeline steps.'''

COMBINED_STEP_ORDER = [ "rampage-idr" ]
'''The (artifically) linear order of all pipeline steps.'''

REP_STEPS = {
    "concatR1": {
                "app":     "concat-fastqs",
                "params":  { "concat_id":  "concat_id" },
                "inputs":  { "reads1_set": "reads_set" },
                "results": { "reads1":     "reads"     }
    },
    "concatR2": {
                "app":     "concat-fastqs",
                "params":  { "concat_id2": "concat_id" },
                "inputs":  { "reads2_set": "reads_set" },
                "results": { "reads2":     "reads"     }
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
                     "rampage_peaks_bb":  "rampage_peaks_bb",
                     "rampage_peaks_gff": "rampage_peaks_gff" }
    }
}
COMBINED_STEPS = {
    "rampage-idr": {
        "app":     "rampage-idr",
        "params":  {},
        "inputs":  { "peaks_a": "peaks_a", 
                     "peaks_b": "peaks_b", 
                     "chrom_sizes": "chrom_sizes" },
        "results": { "rampage_idr_bed": "rampage_idr_bed",
                     "rampage_idr_bb":  "rampage_idr_bb" }
    }
}

FILE_GLOBS = {
    "reads1":               "/*_reads_concat.fq.gz",
    "reads2":               "/*_reads2_concat.fq.gz",
    "all_plus_bw":        "/*_rampage_5p_plusAll.bw",
    "rampage_marked_bam": "/*_rampage_star_marked.bam",
    "all_minus_bg":       "/*_rampage_5p_minusAll.bg",
    "unique_plus_bg":     "/*_rampage_5p_plusUniq.bg",
    "rampage_peaks_bed":  "/*_rampage_peaks.bed",
    "rampage_peaks_bb":   "/*_rampage_peaks.bb",
    "rampage_peaks_gff":  "/*_rampage_peaks.gff",
    "peaks_a":            "/*_rampage_peaks.gff",
    "peaks_b":            "/*_rampage_peaks.gff",
    "all_plus_bg":        "/*_rampage_5p_plusAll.bg",
    "rampage_star_log":   "/*_rampage_star_Log.final.out",
    "all_minus_bw":       "/*_rampage_5p_minusAll.bw",
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
                    help="Biological replicate number",
                    type=int,
                    default=0,
                    required=False)

    ap.add_argument('--tr', '--technical-replicate',
                    help="Technical replicate number (default: 1)",
                    type=int,
                    default=1,
                    required=False)

    ap.add_argument('-c', '--control',
                    help='The control bam for peak calling.',
                    required=False)

    ap.add_argument('--cr','--combine-replicates',
                    help="Combine or compare two replicates (e.g.'1 2_2').'",
                    nargs='+',
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
    '''Adds common and pipeline specific variables to a dict, for use building the workflow.'''
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
    psv = dxencode.common_variables(args,RESULT_FOLDER_DEFAULT,controls=True)
    if psv['exp_type'] != 'rampage':
        print "Experiment %s is not for rampage but for '%s'" % (psv['experiment'],psv['exp_type'])
        sys.exit(1)
    
    # Could be multiple annotations supported per genome
    psv['annotation'] = args.annotation
    if psv['genome'] != GENOME_DEFAULT and psv['annotation'] == ANNO_DEFAULT:
        psv['annotation'] = ANNO_DEFAULTS[psv['genome']]
    if psv['annotation'] not in ANNO_ALLOWED[psv['genome']]:
        print psv['genome']+" has no "+psv['annotation']+" annotation."
        sys.exit(1)

    if not psv['paired_end']:
        print "Rampage is always expected to be paired-end but mapping says otherwise."
        print mapping
        sys.exit(1)

    # Some specific settings
    psv['nthreads']   = 8
    
    # run will either be for combined or single rep.
    if not psv['combined']:
        run = psv['reps']['a']  # If not combined then run will be for the first (only) replicate
    else:
        run = psv
        
    # workflow labeling
    psv['description'] = "The ENCODE Rampage RNA pipeline for long RNAs"
    run['name'] = "rampage_"+psv['genome']
    if psv['genome'] == 'mm10':
        run['name'] += psv['annotation']
    if psv['gender'] == 'female':
        run['name'] += "XX"
    else:
        run['name'] += "XY"
    run['title'] = "Rampage RNA " + psv['experiment'] + " - " + run['rep_tech']
    run['name'] += "_"+psv['experiment']+"_" + run['rep_tech']
    if not psv['combined']:
        run['title'] += " [library '"+run['library_id']+"']"
    run['title'] += " on " + psv['genome']+" - "+psv['gender']

    if verbose:
        print "Pipeline Specific Vars:"
        print json.dumps(psv,indent=4)
    return psv


def find_ref_files(priors,psv):
    '''Locates all reference files based upon organism and gender.'''
    refFiles = {}
    starIx = psv['refLoc']+GENOME_REFERENCES['star_index'][psv['genome']][psv['gender']][psv['annotation']]
    starIxFid = dxencode.find_file(starIx,dxencode.REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    anno = psv['refLoc']+GENOME_REFERENCES['gene_annotation'][psv['genome']][psv['annotation']]
    anno_fid = dxencode.find_file(anno,dxencode.REF_PROJECT_DEFAULT)
    if anno_fid == None:
        sys.exit("ERROR: Unable to locate Gene Annotation file '" + anno + "'")
    else:
        priors['gene_annotation'] = anno_fid

    chromSizes = psv['refLoc']+GENOME_REFERENCES['chrom_sizes'][psv['genome']][psv['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,dxencode.REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chrom_sizes'] = chromSizesFid
    psv['ref_files'] = GENOME_REFERENCES.keys()
    

def find_control_file(psv,rep,project,default=None):
    '''Attempts to find an appropriate control file.'''
    # TODO Make more generic and move to dxencode.py when needed.
    
    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    for file_key in rep['controls']:
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
    print "Unable to find control in search of %s" % rep['controls']
    sys.exit(1)
            

def find_combined_inputs(steps,psv,file_globs,proj_id):
    '''Finds the inputs for a combined run in the directories of the replicate runs.'''
    # TODO Make more generic and move to dxencode.py when needed.
    
    inputs = {}
    for step in psv['path']:
        for fileToken in steps[step]['inputs'].keys():
            rep_key = fileToken[-1]
            if rep_key not in psv['reps']:
                continue
            rep = psv['reps'][rep_key]
            # TODO: No need for multiples at this time.  Deal with it when it comes up.
            fid = dxencode.find_file(rep['resultsFolder'] + file_globs[fileToken],proj_id, \
                                                                      multiple=False, recurse=False)
            if fid != None:
                psv['priors'][fileToken] = fid
                inputs[fileToken] = [ fid ]
            else:
                print "Error: Necessary '%s' for combined run, not found in '%s'." \
                                                                % (fileToken, rep['resultsFolder'])
                print "       Please run for single replicate first."
                sys.exit(1)
    return inputs
    

#######################
def main():

    args = get_args()
    psv = pipeline_specific_vars(args)

    project = dxencode.get_project(psv['project'])
    projectId = project.get_id()

    #print "Building apps dictionary..."
    #pipeSteps, file_globs = dxencode.build_simple_steps(REP_STEP_ORDER,projectId,verbose=True)
    pipeRepSteps = REP_STEPS
    file_globs = FILE_GLOBS
    for rep in psv['reps'].values():
        rep['path'] = REP_STEP_ORDER
    pipeCombinedSteps = None
    if psv['combined']:
        psv['path'] = COMBINED_STEP_ORDER
        pipeCombinedSteps, combined_globs = dxencode.build_simple_steps(psv['path'],projectId)
        file_globs.update(combined_globs)

    # finding fastqs and prior results in a stadardized way
    dxencode.finding_rep_inputs_and_priors(psv,pipeRepSteps,file_globs,project,args.test)

    print "Checking for control files..."
    for rep in psv['reps'].values():
        control = find_control_file(psv,rep,project,args.control)
        rep['inputs']['Control'] = dxencode.find_and_copy_read_files(rep['priors'], [ control ], \
                                                args.test, 'control_bam', rep['resultsFolder'], \
                                                False, projectId)

    if psv['combined']:
        print "Checking for combined-replicate inputs..."
        psv['priors'] = {}
        psv['inputs'] = find_combined_inputs(pipeCombinedSteps,psv,file_globs,projectId)

    # finding pipeline specific reference files in a stadardized way
    dxencode.find_all_ref_files(psv,find_ref_files)

    # deterine steps to run in a stadardized way
    dxencode.determine_steps_needed(psv, pipeRepSteps, pipeCombinedSteps, projectId, args.force)

    # Preperation is done. At this point on we either run rep 'a' or combined.
    if not psv['combined']:
        run = psv['reps']['a']
        run['steps'] = pipeRepSteps
    else:
        run = psv
        run['steps'] = pipeCombinedSteps
        
    dxencode.report_build_launch(psv, run, projectId, test=args.test, launch=args.run)
            
    print "(success)"

if __name__ == '__main__':
    main()

