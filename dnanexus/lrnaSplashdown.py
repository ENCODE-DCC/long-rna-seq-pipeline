#!/usr/bin/env python
import argparse
import os
import sys
import json
import itertools

import dxpy
import dxencode as dxencode
#from dxencode import dxencode as dxencode
import json

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

GENOMES_SUPPORTED = ['hg19', 'mm10']
GENOME_DEFAULT = 'hg19'
''' This the default Genome that long RNA-seq experiments are mapped to.'''

ANNO_DEFAULTS = {'hg19': 'v19', 'mm10': 'M4' }
ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                 'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
ANNO_DEFAULT = ANNO_DEFAULTS[GENOME_DEFAULT]
''' Multiple annotations might be supported for each genome.'''

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/'
''' This the default folder that reference files are found in.'''

RESULT_FOLDER_DEFAULT = '/lrna/'
''' This the default location to place results folders for each experiment.'''

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

REP_STEP_ORDER = {
    # for SE or PE the list in order of steps to run
    "se": [ "concatR1",             "align-tophat-se", "topBwSe", "align-star-se", "starBwSe", "quant-rsem" ],
    "pe": [ "concatR1", "concatR2", "align-tophat-pe", "topBwPe", "align-star-pe", "starBwPe", "quant-rsem" ]
    # examples for testing build_simple_steps
    #"se": [ "align-tophat-se", "align-star-se", "quant-rsem" ],
    #"pe": [ "align-tophat-pe", "align-star-pe", "quant-rsem" ]
    }
'''The (artifically) linear order of all pipeline steps for single or paired-end.'''

REP_STEPS = {
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
        "output_type": "unique minus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_plus_all_bw":   {
        "file_format": "bigWig",
        "output_type": "multi-read plus signal",
        "derived_from": ["tophat_bam"]
    },
    "tophat_plus_uniq_bw":  {
        "file_format": "bigWig",
        "output_type": "unique minus signal",
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
    "rsem_gene_results":    {
        "file_format": "tsv",
        "output_type": "genome quantifications",
        "derived_from": []
        # note should be derived from star_anno_bam
    }
}

TNX_TEMPLATES = {
    "star_anno_bam":        {
        "file_format": "bam",
        "output_type": "transcriptome alignments",
        "derived_from": ["reads1", "reads2"]
    },
    "rsem_iso_results":     {
        "file_format": "tsv",
        "output_type": "transcript quantifications",
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
                                            REF_PROJECT_DEFAULT + ":" + REF_FOLDER_DEFAULT + "')",
                    default=REF_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--testserver',
                    help="Use the test server designated in keypairs.json",
                    action='store_true',
                    required=False)

    ap.add_argument('--test',
                    help='Test run only, do not launch anything.',
                    action='store_true',
                    required=False)

    ap.add_argument('--skipvalidate',
                    help='Skip running Validate Files',
                    action='store_true',
                    required=False)

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
    
    # Could be multiple annotations supported per genome
    psv['annotation'] = args.annotation
    if psv['genome'] != GENOME_DEFAULT and psv['annotation'] == ANNO_DEFAULT:
        psv['annotation'] = ANNO_DEFAULTS[psv['genome']]
    if psv['annotation'] not in ANNO_ALLOWED[psv['genome']]:
        print psv['genome']+" has no "+psv['annotation']+" annotation."
        sys.exit(1)
    
    # Some specific settings
    psv['nthreads']   = 8
    psv['rnd_seed']   = 12345

    # run will either be for combined or single rep.
    if not psv['combined']:
        run = psv['reps']['a']  # If not combined then run will be for the first (only) replicate
    else:
        run = psv
        print "Long-RNA-seq pipeline currently does not support combined-replicate processing."
        print mapping
        sys.exit(1)

    # workflow labeling
    psv['description'] = "The ENCODE RNA Seq pipeline for long RNAs"
    run['name'] = "lrna_"+psv['genome']
    if psv['genome'] == 'mm10':
        run['name'] += psv['annotation']
    if psv['gender'] == 'female':
        run['name'] += "XX"
    else:
        run['name'] += "XY"
    if psv['paired_end']:
        run['title'] = "long RNA-seq paired-end "
        run['name'] += "PE"
    else:
        run['title'] = "long RNA-seq single-end "
        run['name'] += "SE"
    run['title']   += psv['experiment']+" - "+run['rep_tech'] + " (library '"+run['library_id']+"')"
    run['subTitle'] = psv['genome']+", "+psv['gender']+" and annotation '"+psv['annotation']+"'."
    run['name']    += "_"+psv['experiment']+"_"+run['rep_tech']

    if verbose:
        print "Pipeline Specific Vars:"
        print json.dumps(psv,indent=4)
    return psv


def find_ref_files(priors,psv):
    '''Locates all reference files based upon gender, organism and annotation.'''
    refFiles = {}
    topIx = psv['refLoc']+GENOME_REFERENCES['tophat_index'][psv['genome']][psv['gender']][psv['annotation']]
    topIxFid = dxencode.find_file(topIx,dxencode.REF_PROJECT_DEFAULT)
    if topIxFid == None:
        sys.exit("ERROR: Unable to locate TopHat index file '" + topIx + "'")
    else:
        priors['tophat_index'] = topIxFid

    starIx = psv['refLoc']+GENOME_REFERENCES['star_index'][psv['genome']][psv['gender']][psv['annotation']]
    starIxFid = dxencode.find_file(starIx,dxencode.REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['star_index'] = starIxFid

    rsemIx = psv['refLoc']+GENOME_REFERENCES['rsem_index'][psv['genome']][psv['annotation']]
    rsemIxFid = dxencode.find_file(rsemIx,dxencode.REF_PROJECT_DEFAULT)
    if rsemIxFid == None:
        sys.exit("ERROR: Unable to locate RSEM index file '" + rsemIx + "'")
    else:
        priors['rsem_index'] = rsemIxFid

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

    project = dxencode.get_project(psv['project'])
    projectId = project.get_id()

    #print "Building apps dictionary..."
    pipeRepPath = REP_STEP_ORDER['se']
    if psv['paired_end']:
        pipeRepPath = REP_STEP_ORDER['pe']
    #pipeRepSteps, file_globs = dxencode.build_simple_steps(pipeRepPath,projectId,verbose=True)
    pipeRepSteps = REP_STEPS
    file_globs = FILE_GLOBS
    for rep in psv['reps'].values():
        rep['path'] = pipeRepPath

    # finding fastqs and prior results in a stadardized way
    dxencode.finding_rep_inputs_and_priors(psv,pipeRepSteps,file_globs,project,args.test)

    # finding pipeline specific reference files in a stadardized way
    dxencode.find_all_ref_files(psv,find_ref_files)

    # deterine steps to run in a stadardized way
    dxencode.determine_steps_needed(psv, pipeRepSteps, None, projectId)

    #run = psv['reps']['a']
    #run['steps'] = pipeRepSteps
    #dxencode.report_run_plans(psv, run)

    for rep in psv['reps'].values():
        if rep['paired_end']: 
            rep['accessions'] = { 'reads1': [], 'reads2': [] }
        else:
            rep['accessions'] = { 'reads1': [] }
        for rn in rep['fastqs'].keys():
            token = 'reads' + rn
            for fn in rep['fastqs'][rn]:
                rep['accessions'][token].append(fn.split('.')[0]) # split the accession off filename
    if psv['combined']:
        psv['accessions'] = { 'reads1': [], 'reads2': [] }
                
    exp = dxencode.get_exp(psv['experiment'])

    # TODO: Prevent resubmissions
    for run in  (psv['reps'].values() + [ psv ]):
        if 'combined' in run and not psv['combined']:
            continue # only do combined if it was set up

        if len(run['stepsToDo']):
            print "Pipeline %s:%s incomplete, please resubmit jobs: %s" % \
                                (psv['experiment'],run['rep_tech'],run['stepsToDo'])
    
        print "Checking for currently running analyses..."
        dxencode.check_run_log(rep['resultsFolder'],projectId, verbose=True)

        to_submit = [ k for k in run['priors'].keys() if POST_TEMPLATES.get(k) ]
        n = 0 # skip reads
        print "Attempting to submit %s files to args.experiment" % len(to_submit)
        while(to_submit):
            if n > len(run['priors']) * len(run['priors']):
                print "Too many itereations: %s" % run['priors']
                break
            token = to_submit.pop(0)
            print "%s %s" % (token, run['priors'][token])
            f_ob = POST_TEMPLATES.get(token, None)
            n += 1
            if f_ob:
                derive_check = f_ob.get('derived_from', [])
                if derive_check:
                    derived = [ run['accessions'][f] for f in derive_check if run['accessions'].get(f) ]
                    if not derived:
                        to_submit.append(token)
                        continue
                    else:
                        f_ob['derived_from'] = list(itertools.chain(*derived))
                dxFile = dxpy.DXFile(dxid=run['priors'][token])
                print "Post File: %s %s" % (token, dxFile.name)
                f_ob['dataset'] = args.experiment
                f_ob['lab'] = exp['lab']['@id']
                f_ob['award'] = exp['award']['@id']
                f_ob['assembly'] = psv['genome']
                f_ob['genome_annotation'] = psv['annotation']
                ## temporary haxors until file display works
                if 'replicate_id' not in run:
                    f_ob['replicate'] = rep['replicate_id']
                f_ob['notes'] = json.dumps(dxencode.create_notes(dxFile, dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)')))
                print json.dumps(f_ob, sort_keys=True, indent=4, separators=(',',': '))
                if args.test:
                    fake_acc = 'ENCFF%03dAAA' % n
                    print "Fake submission: %s" % fake_acc
                    run['accessions'][token] = [ fake_acc ]
                else:
                    f_ob['derived_from'] = list(itertools.chain(*derived))
            dxFile = dxpy.DXFile(dxid=priors[token])
            print "Post File: %s %s" % (token, dxFile.name)
            f_ob['dataset'] = args.experiment
            f_ob['lab'] = '/labs/j-michael-cherry/'
            f_ob['award'] = '/awards/U41HG006992/'
            f_ob['assembly'] = mapping['genome']
            f_ob['genome_annotation'] = args.annotation
            ## temporary haxors until file display works
            f_ob['replicate'] = mapping['replicate_id']
            f_ob['notes'] = json.dumps(dxencode.create_notes(dxFile, dxencode.get_sw_from_log(dxFile, '\* (\S+)\s+version:\s+(\S+)')))
            print json.dumps(f_ob, sort_keys=True, indent=4, separators=(',',': '))
            if args.testserver:
                server = 'test'
            else:
                server = 'www'

            if args.test:
                fake_acc = 'ENCFF%03dAAA' % n
                print "Fake submission: %s" % fake_acc
                submitted[token] = [ fake_acc ]
            else:
                applet = dxencode.find_applet_by_name('validate-post', projectId )
                job = applet.run({
                    "pipe_file": dxpy.dxlink(dxFile),
                    "file_meta": f_ob,
                    "key": server,
                    "debug": True,
                    "skipvalidate": args.skipvalidate or False
                    })
                print "Submitting %s" % job.id
                job.wait_on_done(interval=1)
                accession = job.describe()['output'].get('accession', "Unknown Acc")
                error = job.describe()['output'].get('error', "Unknown Error")
                submitted[token] = [ accession ]
                print "Posted (%s): %s" % (error, accession)

    # Exit if test only
    if args.test:
        print "Fake submitted %s files." % n
    if args.test:
        sys.exit(0)


if __name__ == '__main__':
    main()

