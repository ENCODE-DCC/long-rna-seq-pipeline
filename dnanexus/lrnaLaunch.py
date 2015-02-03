#!/usr/bin/env python
# lrnaLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
#from dxencode import dxencode as dxencode
import dxencode as dxencode
from launch import Launch

class LrnaLaunch(Launch):
    '''Descendent from Launch class with 'long-rna-seq' methods'''

    PIPELINE_NAME = "long-rna-seq"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline analysis for one replicate. "
    ''' This pipline does not support combined replicates.'''
                    
    ANNO_DEFAULTS = {'hg19': 'v19', 'mm10': 'M4' }
    ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                     'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
    ANNO_DEFAULT = ANNO_DEFAULTS[Launch.GENOME_DEFAULT]
    ''' Multiple annotations might be supported for each genome.'''

    RESULT_FOLDER_DEFAULT = '/lrna/'
    ''' This the default location to place results folders for each experiment.'''

    REP_STEP_ORDER = {
        "se": [ "concatR1",             "align-tophat-se", "topBwSe", "align-star-se", "starBwSe", "quant-rsem" ],
        "pe": [ "concatR1", "concatR2", "align-tophat-pe", "topBwPe", "align-star-pe", "starBwPe", "quant-rsem" ]
        }
    '''The (artifically) linear order of all pipeline steps for single or paired-end.'''
    COMBINED_STEP_ORDER = []
    '''No combined steps in this pipeline.'''

    REP_STEPS = {
        # for each step: app, list of any params, inputs and results
        # TODO: Any results files not listed here would not be 'deprecated' on reruns.
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
                    "params":  { "paired_end":        "paired_end" },  #, "nthreads", "rnd_seed"
                    "inputs":  { "star_anno_bam":     "star_anno_bam",
                                 "rsem_index":        "rsem_index" },
                    "results": { "rsem_iso_results":  "rsem_iso_results",
                                 "rsem_gene_results": "rsem_gene_results" }
                    }
        }

    FILE_GLOBS = {
        # For looking up previous result files, use wild-cards
        "reads1":               "/*_reads_concat.fq.gz",
        "reads2":               "/*_reads2_concat.fq.gz",
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

    REFERENCE_FILES = {
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

    def __init__(self):
        Launch.__init__(self)
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self,parse=False)
        
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
        
        # Some specific settings
        psv['nthreads']   = 8
        psv['rnd_seed']   = 12345

        # run will either be for combined or single rep.
        if not psv['combined']:
            run = psv['reps']['a']  # If not combined then run will be for the first (only) replicate
        else:
            run = psv
            print "Long-RNA-seq pipeline currently does not support combined-replicate processing."
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
        '''Locates all reference files based upon gender, organism and annotation.'''
        topIx = self.psv['refLoc']+self.REFERENCE_FILES['tophat_index'][self.psv['genome']][self.psv['gender']][self.psv['annotation']]
        topIxFid = dxencode.find_file(topIx,dxencode.REF_PROJECT_DEFAULT)
        if topIxFid == None:
            sys.exit("ERROR: Unable to locate TopHat index file '" + topIx + "'")
        else:
            priors['tophat_index'] = topIxFid

        starIx = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][self.psv['gender']][self.psv['annotation']]
        starIxFid = dxencode.find_file(starIx,dxencode.REF_PROJECT_DEFAULT)
        if starIxFid == None:
            sys.exit("ERROR: Unable to locate STAR index file '" + starIx + "'")
        else:
            priors['star_index'] = starIxFid

        rsemIx = self.psv['refLoc']+self.REFERENCE_FILES['rsem_index'][self.psv['genome']][self.psv['annotation']]
        rsemIxFid = dxencode.find_file(rsemIx,dxencode.REF_PROJECT_DEFAULT)
        if rsemIxFid == None:
            sys.exit("ERROR: Unable to locate RSEM index file '" + rsemIx + "'")
        else:
            priors['rsem_index'] = rsemIxFid

        chromSizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chromSizesFid = dxencode.find_file(chromSizes,dxencode.REF_PROJECT_DEFAULT)
        if chromSizesFid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
        else:
            priors['chrom_sizes'] = chromSizesFid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()


    #######################

if __name__ == '__main__':
    '''Run from the command line.'''
    lrnaLaunch = LrnaLaunch()
    lrnaLaunch.run()

