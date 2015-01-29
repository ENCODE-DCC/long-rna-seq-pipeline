#!/usr/bin/env python
# rampageLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
#from dxencode import dxencode as dxencode
import dxencode as dxencode
from launch import Launch

class RampageLaunch(Launch):
    '''Descendent from Launch class with 'rampage' methods'''

    PIPELINE_NAME = "rampage"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''
    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    ANNO_DEFAULTS = {'hg19': 'v19', 'mm10': 'M4' }
    ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                     'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
    ANNO_DEFAULT = ANNO_DEFAULTS[Launch.GENOME_DEFAULT]
    ''' Multiple annotations might be supported for each genome.'''

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

    REFERENCE_FILES = {
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


    def __init__(self):
        Launch.__init__(self)
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self,parse=False)
        
        # NOTE: Could override get_args() to have this non-generic control message
        #ap.add_argument('-c', '--control',
        #                help='The control bam for peak calling.',
        #                required=False)

        ap.add_argument('-a', '--annotation',
                        help="Label of annotation (default: '" + self.ANNO_DEFAULT + "')",
                        choices=[self.ANNO_DEFAULT, 'M2','M3','M4'],
                        default=self.ANNO_DEFAULT,
                        required=False)
        return ap.parse_args()

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self,args)
        
        # Could be multiple annotations supported per genome
        psv['annotation'] = args.annotation
        if psv['genome'] != self.GENOME_DEFAULT and psv['annotation'] == self.ANNO_DEFAULT:
            psv['annotation'] = self.ANNO_DEFAULTS[psv['genome']]
        if psv['annotation'] not in self.ANNO_ALLOWED[psv['genome']]:
            print psv['genome']+" has no "+psv['annotation']+" annotation."
            sys.exit(1)

        if not psv['paired_end']:
            print "Rampage is always expected to be paired-end but mapping says otherwise."
            sys.exit(1)

        # Some specific settings
        psv['nthreads']   = 8
        psv['control'] = args.control
        
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

        # Must override results location because of annotation
        psv['resultsLoc'] = args.resultsLoc
        if psv['resultsLoc'] == self.RESULT_FOLDER_DEFAULT:
            if psv['genome'] == 'mm10' and 'annotation' in psv:
                psv['resultsLoc'] = self.RESULT_FOLDER_DEFAULT + psv['genome'] + '/' + psv['annotation'] + '/'
            else:
                psv['resultsLoc'] = self.RESULT_FOLDER_DEFAULT + psv['genome'] + '/'
        if not psv['resultsLoc'].endswith('/'):
            psv['resultsLoc'] += '/' 
        psv['resultsFolder'] = psv['resultsLoc'] + psv['experiment'] + '/'
        psv['reps']['a']['resultsFolder'] = psv['resultsLoc'] + psv['experiment'] + '/' + \
                                                              psv['reps']['a']['rep_tech'] + '/'
        if psv['combined']:
            psv['reps']['b']['resultsFolder'] = psv['resultsLoc'] + psv['experiment'] + '/' + \
                                                                  psv['reps']['b']['rep_tech'] + '/'

        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        starIx = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][self.psv['gender']][self.psv['annotation']]
        starIxFid = dxencode.find_file(starIx,dxencode.REF_PROJECT_DEFAULT)
        if starIxFid == None:
            sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
        else:
            priors['star_index'] = starIxFid

        anno = self.psv['refLoc']+self.REFERENCE_FILES['gene_annotation'][self.psv['genome']][self.psv['annotation']]
        anno_fid = dxencode.find_file(anno,dxencode.REF_PROJECT_DEFAULT)
        if anno_fid == None:
            sys.exit("ERROR: Unable to locate Gene Annotation file '" + anno + "'")
        else:
            priors['gene_annotation'] = anno_fid

        chromSizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chromSizesFid = dxencode.find_file(chromSizes,dxencode.REF_PROJECT_DEFAULT)
        if chromSizesFid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
        else:
            priors['chrom_sizes'] = chromSizesFid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
    

    def find_control_file(self,rep,default=None):
        '''Attempts to find an appropriate control file.'''
        # TODO Make more generic and move to dxencode.py when needed.
        
        (AUTHID,AUTHPW,SERVER) = dxencode.processkey(self.server_key)
        for file_key in rep['controls']:
            url = '%s%s/?format=json&frame=embedded' % (SERVER,file_key)
            #print '-- ' + AUTHID + " " + AUTHPW + " " + SERVER + " " + url
            try:
                response = dxencode.encoded_get(url, AUTHID, AUTHPW)
                file_obj = response.json()
            except:
                print "URL to control [%s] returned ?" % url
                print response
                sys.exit(1)
            #print json.dumps(response,indent=4)
            rep_id = file_obj["replicate"]['@id']
            url = '%s%s/?format=json&frame=embedded' % (SERVER,rep_id)
            try:
                response = dxencode.encoded_get(url, AUTHID, AUTHPW)
                rep_obj = response.json()
            except:
                print "URL to replicate [%s] returned ?" % url
                print response
                sys.exit(1)
            exp_id = rep_obj['experiment'].split('/')[2]
            rep_tech = "rep%s_%s" % \
                    (rep_obj['biological_replicate_number'], rep_obj['technical_replicate_number'])
            # default by cheating
            control_root = self.CONTROL_ROOT_FOLDER
            path_n_glob = control_root + exp_id + '/' + rep_tech + '/' + self.CONTROL_FILE_GLOB
            target_folder = dxencode.find_folder(exp_id + '/' + rep_tech,self.project,control_root)
            #print "Target found [%s]" % target_folder
            if target_folder != None:
                path_n_glob = target_folder + '/' + self.CONTROL_FILE_GLOB
            fid = dxencode.find_file(path_n_glob,self.proj_id,multiple=False,recurse=False)
            if fid != None:
                return dxencode.file_path_from_fid(fid)
                
        if default != None:
            return default
        print "Unable to find control in search of %s" % rep['controls']
        sys.exit(1)
            

    #######################

if __name__ == '__main__':
    '''Run from the command line.'''
    rampageLaunch = RampageLaunch()
    rampageLaunch.run()

