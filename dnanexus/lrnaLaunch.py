#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
from datetime import datetime

import dxpy

# NOTES: This command-line utility will run the long RNA-seq pipeline for a single replicate
#      - All results will bw written to a folder /runs/<expId>/rep<#>.
#      - If any results are already in that directory, then the steps that created those results 
#        will not be rerun.  
#      - If any jobs for the experiment and replicate are already running, nothing new will be 
#        launched.
#      - Most (BUT NOT ALL) of the code is generic, relying on the hard-coded JSON below and
#        relying on hard-coded tokens for step, file and parameter names:
#        - STEP_ORDER is the list of steps for single or paired-end reads
#        - STEPS contains step definitions and enforces dependencies by result file tokens matching
#          input file tokens of later steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.
#        Long RNA-seq specific code is marked with '### LRNA specific'

GENOME_DEFAULT = 'hg19'
''' This the default Genome that long RNA-seq experiments are mapped to.'''

ANNO_DEFAULT = 'v19'

PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'scratchPad'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/ref'
''' This the default folder that reference files are found in.'''

RESULT_FOLDER_DEFAULT = '/runs'
''' This the default location to place results folders for each experiment.'''

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

STEP_ORDER = {
    # for SE or PE the list in order of steps to run
    'se': [ 'concatR1',             'tophatSe', 'topBwSe', 'starSe', 'starBwSe', 'rsem' ],
    'pe': [ 'concatR1', 'concatR2', 'tophatPe', 'topBwPe', 'starPe', 'starBwPe', 'rsem' ]
    }

STEPS = {
    # for each step: app, list of any params, inputs and results (both as fileToken: app_obj_name)
    'concatR1': {
                'app':     'concat-fastqs',
                'params':  { 'rootR1': 'outfile_root' },
                'inputs':  { 'reads1_set': 'fastq_files' },
                'results': { 'reads1': 'combined_fastq' }
                },
    'concatR2': {
                'app':     'concat-fastqs',
                'params':  { 'rootR2': 'outfile_root' },
                'inputs':  { 'reads2_set': 'fastq_files' },
                'results': { 'reads2': 'combined_fastq' }
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
    'rsemIso':         "/*isoforms.results",
    'rsemGene':        "/*genes.results"
    }

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should remove annotation if only one per genome
    # TODO: should use ACCESSION based fileNames
    'tophatIndex':  {
                    'hg19': {
                            'female':   {
                                        'v19':          'female_hg19_v19_combined_tophatIndex.tgz',
                                        'v19_annoOnly': 'female_hg19_v19_annoOnly_tophatIndex.tgz'
                                        },
                            'male':     {
                                        'v19':          'male_hg19_v19_combined_tophatIndex.tgz',
                                        'v19_annoOnly': 'male_hg19_v19_annoOnly_tophatIndex.tgz'
                                        }
                            },
                    'mm10': {
                            'female':   { 'M3':         'mm10_female_M3_tophatIndex.tgz' },
                            'male':     { 'M3':         'mm10_male_M3_tophatIndex.tgz'   }
                            }
                    },
    'starIndex':    {
                    'hg19': {
                            'female':   {
                                        'v19':          'female_hg19_v19_combined_starIndex.tgz',
                                        'v19_annoOnly': 'female_hg19_v19_annoOnly_starIndex.tgz'
                                        },
                            'male':     {
                                        'v19':          'male_hg19_v19_combined_starIndex.tgz',
                                        'v19_annoOnly': 'male_hg19_v19_annoOnly_starIndex.tgz'
                                        }
                            },
                    'mm10': {
                            'female':   { 'M3':         'mm10_female_M3_starIndex.tgz' },
                            'male':     { 'M3':         'mm10_male_M3_starIndex.tgz'   }
                            }
                    },
    'rsemIndex':    {
                    'hg19': {
                            'v19':          'male_hg19_v19_combined_rsemIndex.tgz',
                            'v19_annoOnly': 'male_hg19_v19_annoOnly_rsemIndex.tgz'
                            },
                    'mm10': { 'M3':         'mm10_male_M3_rsemIndex.tgz' }
                    },
    'chromSizes':   {
                    'hg19': {
                            'female':   'female_hg19.chrom.sizes',
                            'male':     'male_hg19.chrom.sizes'
                            },
                    'mm10': {
                            'female':   'mm10/mm10_female.chrom.sizes',
                            'male':     'mm10/mm10_male.chrom.sizes'
                            }
                    }
    }

APPLETS = {}
# Used for caching applets that might be called more than once in pipeline
FILES = {}
# Used for caching file dxlinks that might be needed more than once in building the workflow

def get_args():
    '''Parse the input arguments.'''
    ### LRNA specific
    ap = argparse.ArgumentParser(description="Launches long RNA-seq pipeline analysis for " + 
                "one replicate on single or paired-end reads. Can be run repeatedly and will " + 
                "launch only the steps that are needed to finished the pipeline. All results " + 
                "will be placed in the folder /<resultsLoc>/<experiment>/<replicate>.")
    ### LRNA specific

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    required=True)

    ap.add_argument('-r', '--replicate',
                    help="Replicate number (default: 1)",
                    choices=['1', '2'],
                    default='1',
                    required=True)

    ap.add_argument('-1', '--reads1',
                    help='Only reads fastq file or first of pair-end reads.',
                    nargs='+',
                    required=True)

    ap.add_argument('-2', '--reads2',
                    help='The second of paired-end reads files.',
                    nargs='+',
                    required=False)

    ### LRNA specific
    ap.add_argument('-l', '--library',
                    help='ENCODE accession of biosample library (for BAM header)',
                    required=True)
    ### LRNA specific

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

    ### LRNA specific
    # TODO: should remove annotation if only one per genome
    ap.add_argument('-a', '--annotation',
                    help="Label of annotation (default: '" + ANNO_DEFAULT + "')",
                    choices=[ANNO_DEFAULT, 'v19_annoOnly','M3'],
                    default=ANNO_DEFAULT,
                    required=False)
    ### LRNA specific

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    action='store_true',
                    default=PROJECT_DEFAULT,
                    required=False)

    ap.add_argument('--refLoc',
                    help="The location to find reference files (default: '" + \
                                            REF_PROJECT_DEFAULT + ":" + REF_FOLDER_DEFAULT + "')",
                    action='store_true',
                    default=REF_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    action='store_true',
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('-t', '--test',
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

def pipelineSpecificExtras(organism, gender, annotation, experiment, replicate, library, pairedEnd):
    '''Adds pipeline specific variables to a dict, for use building the workflow.'''
    ### LRNA specific
    extras = {}
    extras['organism']   = organism
    extras['gender']     = gender
    extras['annotation'] = annotation
    if organism == 'hg19' and annotation not in [ANNO_DEFAULT, 'v19_annoOnly']:
        print organism+" has no "+annotation+" annotation."
        sys.exit(1)
    if organism == 'mm10' and annotation not in ['M3']:
        print organism+" has no '"+annotation+"' annotation."
        sys.exit(1)
    extras['experiment'] = experiment
    extras['replicate']  = replicate
    extras['library']    = library
    extras['pairedEnd']  = pairedEnd 
    # workflow labeling
    genderToken = "XY"
    if gender == 'female':
        genderToken = "XX"
    extras['description'] = "The ENCODE RNA Seq pipeline for long RNAs"
    extras['title'] = "long RNA-seq single-end "
    extras['name'] = "lrna_"+organism+genderToken+"SE_"
    if pairedEnd:
        extras['title'] = "long RNA-seq paired-end "
        extras['name'] = "lrna_"+organism+genderToken+"PE_"
    extras['title'] += experiment+" - rep"+replicate + " (library '"+library+"')"
    extras['name']  += experiment+"_rep"+replicate
    extras['subTitle'] = organism+", "+gender+" and annotation '"+annotation+"'."

    # Non-file app inputs
    extras['rootR1'] = experiment + 'rep' + replicate + '_concatR1'
    extras['rootR2'] = experiment + 'rep' + replicate + '_concatR2'
    ### LRNA specific
    return extras


def getApp(appName, appsProjectId):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '*'
    if (appName, appsProjectId) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=appName,
                                          project=appsProjectId,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(appName, appsProjectId)] = found
        cached = ''

    #print cached + "Resolved %s to %s" % (appName, APPLETS[(appName, appsProjectId)].get_id())
    return APPLETS[(appName, appsProjectId)]


def findFile(filePath,project=None,verbose=False,multiple=False):
    '''Using a DX style file path, find the file.'''
    proj = project
    path = filePath
    fileName = filePath
    if filePath.find(':') != -1:
        proj, path = filePath.split(':', 1)
    if path.rfind('/') != -1:
        path, fileName = path.rsplit('/', 1)
    else:
        fileName = path
        path = ''
    if proj == None:
        if verbose:
            print "ERROR: Don't know what project to use for '" + path + "'."
        return None
    if proj.find('project-') == 0:
        projId = proj
    else:
        projId = getProject(proj, level='VIEW').get_id()
    mode = 'exact'
    if filePath.find('*') or filePath.find('?'):
        mode = 'glob'
    fileDicts = list(dxpy.find_data_objects(classname='file', folder=path, name=fileName, 
                                            name_mode=mode, project=projId, return_handler=False))

    if fileDicts == None or len(fileDicts) == 0:
        #print "- Found 0 files from '" + proj + ":" + filePath + "'."
        if verbose:
            print "ERROR: Failed to find '" + proj + ":" + filePath + "'."
        return None
    elif len(fileDicts) > 1 or multiple:
        #print "- Found "+str(len(fileDict))+" files from '" + proj + ":" + filePath + "'."
        if not multiple:
            if verbose:
                print "ERROR: Found "+str(len(fileDict))+" files when expecting 1 '" + proj + ":" + filePath + "'."
            return None
        fids = []
        for fileDict in fileDicts:
            FILES[fileDict['id']] = dxpy.dxlink(fileDict)
            fids += [ fileDict['id'] ]
        return fids
    else:
        #print "- FOUND '" + proj + ":" + filePath + "'."
        FILES[fileDicts[0]['id']] = dxpy.dxlink(fileDicts[0])
        return fileDicts[0]['id']

def moveFiles(fids, folder, projectId):
    '''Moves files to supplied folder.  Expected to be in the same project.'''
    for fid in fids:
        fileDict = dxpy.describe(FILES[fid]) # FILES contain dxLinks
        if fileDict['project'] != projectId:
            print "ERROR: Failed to move '" + fileDict['name'] + "' as it is not in '" + \
                                                                                projectId + "'."
            sys.exit(1)
    proj = dxpy.DXProject(projectId)
    if not projectFolderExists(proj, folder):
        proj.new_folder(folder,parents=True)
    proj.move(folder,fids)

def copyFiles(fids, projectId, folder, overwrite=False):
    '''Copies array of dx file dicts to project:/folder, returning new array of dx file dicts.'''
    newFids = []
    for fid in fids:
        fileDict = dxpy.describe(FILES[fid]) # FILES contain dxLinks
        if fileDict['project'] == projectId:
            # cannot copy into the same project!!!
            # so just leave in place and pretend that we did!
            #proj = dxpy.DXProject(projectId)
            #proj.move(folder,[fid])
            newFids += [ fid ]
            continue

        # Check to see if file already exists.
        alreadyThere = findFile(folder+'/'+fileDict['name'],projectId)
        if alreadyThere is None or overwrite:
            # remove what is alreadyThere?
            #if alreadyThere is not None:
            #    proj = dxpy.DXProject(projectId)
            #    proj.remove_objects([alreadyThere])
            dxFile = dxpy.get_handler(FILES[fid])
            newLink = dxpy.dxlink(dxFile.clone(projectId, folder))
        else:
            newLink = FILES(alreadyThere)
        if newLink == None:
            print "ERROR: Failed in copy of '" + fileDict['project'] + ":" + fileDict['name'] + \
                    "' to '" + projectId + ":" + folder + "'."
            sys.exit(1)
        newDict = dxpy.describe(newLink)
        FILES[newDict['id']] = newLink
        newFids += [ newDict['id'] ]

    return newFids

def findFileSet(fileSet,projectId=None):
    '''Find all files in a set, and prints error(s) and exits if any are missing.'''
    fids = []
    if len(fileSet): 
        for oneFile in fileSet:
            fid = findFile(oneFile,projectId,verbose=True)
            if fid != None:
                fids += [ fid ]
        if len(fids) != len(fileSet):
            print "ERROR: expecting " + str(len(fileSet)) + " but only found " + str(len(fids)) + "."
            sys.exit(1) # verbose already gave an error message(s)

    return fids

def findTargetFileSet(fileSet,targetFolder,project=None):
    '''Looks for a set of files in a target destination, returning the set.'''
    proj = project
    path = targetFolder
    if targetFolder.find(':') != -1:
        proj, path = targetFolder.split(':')
    if proj.find('project-') == 0:
        projId = proj
    else:
        projId = getProject(proj, level='VIEW').get_id()
    targetFids = []
    for oneFile in fileSet:
        parts = oneFile.rsplit('/',1)
        parts.reverse()
        fid = findFile(path + "/" + parts[0],projId)
        if fid != None:
            targetFids += [ fid ]
    return targetFids

def findAndCopyReadFiles(priors, readSet, testOnly, readToken, resultsFolder, projectId=None):
    '''Looks for read files and copies them to results folder if necessary.'''
    readFiles = []
    if len(readSet) > 0:
        readFiles = findFileSet(readSet,projectId)
        if not testOnly:
            priorReads = findTargetFileSet(readSet,resultsFolder,projectId)
            if len(priorReads) == len(readFiles):
                readFiles = priorReads
            else: # any files to be moved, move all
                if readToken in priors:
                    del priors[readToken] # make sure to regenerate the combined file if necessary
                readFiles = copyFiles(readFiles, projectId, resultsFolder)
        if len(readFiles) > 1:        # If more than 1 file, the 'set' will need the 'concat' step.
            priors[readToken+'_set'] = readFiles
        else:                         # otherwise, the 1 file is same as 'concat' result.
            priors[readToken] = readFiles[0]
    return readFiles 

def findPriorResults(pairedEnd,resultsFolder,projectId):
    '''Looks for all result files in the results folder.'''
    priors = {}
    steps = []
    if pairedEnd:
        steps = STEP_ORDER['pe']
    else:
        steps = STEP_ORDER['se']
    for step in steps:
        for fileToken in STEPS[step]['results'].keys():
            fid = findFile(resultsFolder + FILE_GLOBS[fileToken],projectId)
            if fid != None: 
                priors[fileToken] = fid
    return priors

def findReferenceFiles(priors,refLoc,extras):
    '''Locates all reference files based upon gender, organism and annotation.'''
    refFiles = {}
    ### LRNA specific
    topIx = refLoc+'/'+GENOME_REFERENCES['tophatIndex'][extras['organism']][extras['gender']][extras['annotation']]
    topIxFid = findFile(topIx,REF_PROJECT_DEFAULT)
    if topIxFid == None:
        sys.exit("ERROR: Unable to locate TopHat index file '" + topIx + "'")
    else:
        priors['tophatIndex'] = topIxFid
   
    starIx = refLoc+'/'+GENOME_REFERENCES['starIndex'][extras['organism']][extras['gender']][extras['annotation']]
    starIxFid = findFile(starIx,REF_PROJECT_DEFAULT)
    if starIxFid == None:
        sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
    else:
        priors['starIndex'] = starIxFid

    rsemIx = refLoc+'/'+GENOME_REFERENCES['rsemIndex'][extras['organism']][extras['annotation']]
    rsemIxFid = findFile(rsemIx,REF_PROJECT_DEFAULT)
    if rsemIxFid == None:
        sys.exit("ERROR: Unable to locate RSEM index file '" + rsemIx + "'")
    else:
        priors['rsemIndex'] = rsemIxFid

    chromSizes = refLoc+'/'+GENOME_REFERENCES['chromSizes'][extras['organism']][extras['gender']]
    chromSizesFid = findFile(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chromSizes'] = chromSizesFid
    ### LRNA specific

def determineStepsToDo(pairedEnd, priors, deprecate, projectId, force=False):
    '''Determine what steps need to be done, base upon prior results.'''
    willCreate = []
    stepsToDo = []
    steps = []
    if pairedEnd:
        steps = STEP_ORDER['pe']
    else:
        steps = STEP_ORDER['se']
    for step in steps:
        # Force will include the first step with all its inputs
        # This should avoid forcing concat if it isn't needed
        # 
        if force:
            inputs = STEPS[step]['inputs'].keys()
            count = 0
            for input in inputs:
                if input in priors:
                    count += 1
            if count == len(inputs):
                stepsToDo += [ step ]
        if step not in stepsToDo:
            results = STEPS[step]['results'].keys()
            for result in results:
                if result not in priors:
                    #print "- Adding step '"+step+"' because prior '"+result+"' was not found."
                    stepsToDo += [ step ]
                    break
        # If results are there but inputs are being recreated, then step must be rerun
        if step not in stepsToDo:
            inputs = STEPS[step]['inputs'].keys()
            for inp in inputs:
                if inp in willCreate:
                    #print "- Adding step '"+step+"' due to prior step dependency."
                    stepsToDo += [ step ]
                    break
        # Any step that is rerun, will cause prior results to be deprecated
        # NOTE: It is necessary to remove from 'priors' so succeeding steps are rerun
        # NOTE: It is also important to move prior results out of target folder to avoid confusion!
        if step in stepsToDo:
            results = STEPS[step]['results'].keys()
            for result in results:
                willCreate += [ result ]
                if result in priors:
                    deprecate += [ priors[result] ]
                    del priors[result] 
                    # if results are in folder, then duplicate files cause a problem!
                    # So add to 'deprecate' to move or remove before launching 

    # Now make sure the steps can be found, and error out if not.
    for step in stepsToDo:
        app = STEPS[step]['app']
        dxApp = dxpy.find_data_objects(classname='file', name=app, name_mode='exact', 
                                                         project=projectId, return_handler=False)
        if dxApp == None:
            print "ERROR: failure to locate app '"+app+"'!"
            sys.exit(1)

    return stepsToDo
 
def checkRunsPreviouslyLaunched(resultsFolder,projectId):
    '''Checks for currently running jobs and will exit if found.'''
    launchFilePath = resultsFolder + '/' + RUNS_LAUNCHED_FILE
    launchFids = findFile(launchFilePath,projectId,multiple=True)
    if launchFids == None:
        print "  No prior jobs launched."
    else:
        # NOTE: Appending to the one file, but just in case handle multiple files.
        for fid in launchFids:
            with dxpy.open_dxfile(fid) as fd:
                for line in fd:
                    #print "Looking for job ["+line+"]"
                    runId = line.split(None,1)
                    if not runId[0].startswith('analysis-'):
                        continue
                    analysis = dxpy.DXAnalysis(dxid=runId[0])
                    if analysis == None:
                        continue
                    state = analysis.describe()['state']
                    # states I have seen: in_progress, terminated, done, failed
                    if state not in [ "done", "failed", "terminated" ]:
                        msg="Exiting: Can't launch because prior run ["+runId[0]+"] "
                        if len(runId) > 1:
                            msg+="("+runId[1]+") "
                        msg+= "has not finished (currently '"+state+"')."
                        print msg
                        sys.exit(1)
                    else:
                        msg="  Prior run ["+runId[0]+"] "
                        if len(runId) > 1:
                            msg+="("+runId[1]+") "
                        msg+= "is '"+state+"'."
                        print msg

def logThisRun(runId,resultsFolder,projectId):
    '''Adds a runId to the runsLaunched file in resultsFolder.'''
    # NOTE: DX manual lies?!  Append not possible?!  Then write new/delete old
    launchFilePath = resultsFolder + '/' + RUNS_LAUNCHED_FILE
    oldFid = findFile(launchFilePath,projectId)
    newFh = dxpy.new_dxfile('a',project=projectId,folder=resultsFolder,name=RUNS_LAUNCHED_FILE)
    newFh.write(runId+' started:'+str(datetime.now())+'\n')
    if oldFid is not None:
        with dxpy.open_dxfile(oldFid) as oldFh:
            for oldRunId in oldFh:
                newFh.write(oldRunId+'\n')
        proj = dxpy.DXProject(projectId)
        proj.remove_objects([oldFid])
    newFh.close()

def filePathFromFid(fid,projectToo=False):
    '''Returns full dx path to file from a file id.'''
    fileDict = dxpy.describe(FILES[fid]) # FILES contain dxLinks
    path = fileDict['folder'] + '/' + fileDict['name']
    if projectToo:
        projDict = dxpy.describe(fileDict['project'])
        path = projDict['name'] + ':' + path
    return path

def projectFolderExists(project, folder):
    '''Returns True if folder found in project.'''
    folders = project.list_folder()['folders']
    return folder in folders

def getProject(projectName, level=None):
    '''Returns the DXProject by name or errors out if not found.'''
    try:
        project = dxpy.find_one_project(name=projectName, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print "Could not find 1 and only 1 project named '"+projectName+"'."
        sys.exit(1)

    return dxpy.DXProject(project['id'])


def createWorkflow(stepsToDo, priors, extras, resultsFolder, projectId, appProjectId=None):
    '''This function will populate a workflow for the stepsToDo.'''

    if len(stepsToDo) < 1:
        return None
    if appProjectId == None:
        appProjectId = projectId

    # create a workflow object
    wf = dxpy.new_dxworkflow(title=extras['name'],name=extras['name'],folder=resultsFolder, 
                                            project=projectId,description=extras['description'])

    # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
    prevStepResults = {}
    for step in stepsToDo:
        appName = STEPS[step]['app']
        app = getApp(appName, appProjectId)
        appInputs = {}
        # file inputs
        for fileToken in STEPS[step]['inputs'].keys():
            appInp = STEPS[step]['inputs'][fileToken]
            if fileToken in prevStepResults:
                appInputs[ appInp ] = prevStepResults[fileToken]
            elif fileToken in priors:
                if isinstance(priors[fileToken], list):
                    appInputs[ appInp ] = []
                    for fid in priors[fileToken]:
                        appInputs[ appInp ] += [ FILES[fid] ]
                else:
                    appInputs[ appInp ] = FILES[ priors[fileToken] ]
            else:
                print "ERROR: step '"+step+"' can't find input '"+fileToken+"'!"
                sys.exit(1)
        # Non-file app inputs
        if 'params' in STEPS[step]:
            for param in STEPS[step]['params'].keys():
                appParam = STEPS[step]['params'][param]
                if param in extras:
                    appInputs[ appParam ] = extras[param]
                else:
                    print "ERROR: unable to locate '"+param+"' in extras."
                    sys.exit(1)
        # Add wf stage
        stageId = wf.add_stage(app, stage_input=appInputs, folder=resultsFolder)
        # outputs, which we will need to link to 
        for fileToken in STEPS[step]['results'].keys():
            appOut = STEPS[step]['results'][fileToken]
            prevStepResults[ fileToken ] = dxpy.dxlink({ 'stage': stageId,'outputField': appOut })
    wfRun = wf.run({})
    return wfRun.describe()

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
    extras = pipelineSpecificExtras(args.organism,args.gender,args.annotation, 
                                    args.experiment, args.replicate, args.library, pairedEnd)
    project = getProject(args.project)
    projectId = project.get_id()

    if args.resultsLoc == RESULT_FOLDER_DEFAULT and args.organism == 'mm10':
        args.resultsLoc = RESULT_FOLDER_DEFAULT + '/' + args.organism
    resultsFolder = args.resultsLoc + '/' + args.experiment + '/rep' + args.replicate
    if not args.test:
        if not projectFolderExists(project, resultsFolder):
            project.new_folder(resultsFolder,parents=True)

    print "Checking for prior results..."
    # Check if there are previous results
    # Perhaps reads files are already there?
    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    priors = findPriorResults(pairedEnd,resultsFolder,projectId)

    print "Checking for read files..."
    # Find all reads files and move into place
    reads1 = findAndCopyReadFiles(priors, args.reads1, args.test, 'reads1', resultsFolder, projectId)
    reads2 = findAndCopyReadFiles(priors, args.reads2, args.test, 'reads2', resultsFolder, projectId)

    print "Looking for reference files..."
    findReferenceFiles(priors,args.refLoc,extras)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = determineStepsToDo(pairedEnd, priors, deprecateFiles, projectId, force=args.force)

    # Report the plans
    print "Running '"+extras['title']+"'"
    print "     on "+extras['subTitle']
    if pairedEnd:
        print "- Reads1: "
    else:
        print "- Reads: "
    for fid in reads1:
        print "  " + filePathFromFid(fid)
    if pairedEnd:
        print "- Reads2: "
        for fid in reads2:
            print "  " + filePathFromFid(fid)
    print "- Reference files:"
    for token in GENOME_REFERENCES.keys():
        print "  " + filePathFromFid(priors[token],True)
    print "- Results written to: " + args.project + ":" +resultsFolder
    if len(stepsToDo) == 0:
        print "* All expected results are in the results folder, so there is nothing to do."
        print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
        print "  rerun all steps; or remove suspect results from the folder before launching."
        sys.exit(0)
    else:
        print "- Steps to run:"
        steps = []
        if pairedEnd:
            steps = STEP_ORDER['pe']
        else:
            steps = STEP_ORDER['se']
        for step in steps:
            if step in stepsToDo:
                print "  * "+STEPS[step]['app']+" will be run"
            else:
                if not step.find('concat') == 0:
                    print "    "+STEPS[step]['app']+" has already been run"

    print "Checking for currently running analyses..."
    checkRunsPreviouslyLaunched(resultsFolder,projectId)

    if len(deprecateFiles) > 0:
        if args.test:
            print "Would move "+str(len(deprecateFiles))+" prior result file(s) to '" + \
                                                                    resultsFolder+"/deprecated'."
            for fid in deprecateFiles:
                print "  " + filePathFromFid(fid)
        else:
            print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '" + \
                                                                resultsFolder+"/deprecated'..."
            moveFiles(deprecateFiles,resultsFolder+"/deprecated",projectId)

    # Exit if test only
    if args.test:
        #updateRunsLaunched(extras['experiment'],resultsFolder,projectId)
        print "TEST ONLY - exiting."
        sys.exit(0)

    print "Launch sequence initiating..."
    wfRun = createWorkflow(stepsToDo, priors, extras, resultsFolder,projectId)

    print "  We have liftoff!"
    logThisRun(wfRun['id'],resultsFolder,projectId)

    print "  Launched " + wfRun['id']
    print "(success)"

if __name__ == '__main__':
    main()

