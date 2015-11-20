#!/usr/bin/env python
# srnaLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
from launch import Launch
#from template import Launch # (does not use dxencode at all)

class SrnaLaunch(Launch):
    '''Descendent from Launch class with 'small-rna-seq' methods'''

    PIPELINE_NAME = "small-rna-seq"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline analysis for one replicate. "
    ''' This pipline does not support combined replicates.'''
                    
    ANNO_DEFAULTS = {'hg19': 'v19', 'mm10': 'M4' }
    ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                     'mm10': [ ANNO_DEFAULTS['mm10'] ] }
    ANNO_DEFAULT = ANNO_DEFAULTS[Launch.GENOME_DEFAULT]
    ''' Multiple annotations might be supported for each genome.'''

    PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
    '''Currently there are no combined replicate processing steps.'''
    
    PIPELINE_BRANCHES = {
    #'''Each branch must define the 'steps' and their (artificially) linear order.'''
        "REP": {
                "ORDER": [  "small-rna-align", "small-rna-signals"  ],
                "STEPS": {
                            "small-rna-align": {
                                "app":      "small-rna-align", 
                                "params": { "library_id":   "library_id", 
                                            "nthreads":     "nthreads"      }, 
                                "inputs": { "reads":        "reads", 
                                            "star_index":   "star_index"    }, 
                                "results": {"srna_bam":     "srna_bam",
                                            "srna_quant":   "srna_quant", 
                                            "star_log":     "star_log"      },
                            }, 
                            "small-rna-signals": {
                                "app":      "small-rna-signals", 
                                "inputs": { "srna_bam":         "srna_bam", 
                                            "chrom_sizes":      "chrom_sizes"       }, 
                                "results": {"all_minus_bw":     "all_minus_bw", 
                                            "unique_plus_bw":   "unique_plus_bw", 
                                            "all_plus_bw":      "all_plus_bw", 
                                            "unique_minus_bw":  "unique_minus_bw"   },
                            }, 
                },
        },
        "COMBINED_REPS": {
                "ORDER": [  "small-rna-mad-qc" ],
                "STEPS": {
                            "small-rna-mad-qc": {
                                "app":      "small-rna-mad-qc", 
                                "inputs": { "quants_a":         "quants_a", 
                                            "quants_b":         "quants_b",
                                            "annotations":      "annotations"   }, 
                                "results": {"mad_plot":         "mad_plot"      },
                            },
                
                },
        },
    }
    
    FILE_GLOBS = {
        #"reads":            "/*_concat.fq.gz", 
        "srna_bam":         "/*_srna_star.bam", 
        "star_log":         "/*_srna_star_Log.final.out", 
        "srna_quant":       "/*_srna_star_quant.tsv", 
        "all_minus_bw":     "/*_small_minusAll.bw", 
        "all_plus_bw":      "/*_small_plusAll.bw", 
        "unique_minus_bw":  "/*_small_minusUniq.bw",
        "unique_plus_bw":   "/*_small_plusUniq.bw", 
        "quants_a":         "/*_srna_star_quant.tsv",
        "quants_b":         "/*_srna_star_quant.tsv",
        "mad_plot":         "/*_mad_plot.png", 
    }    

    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should use ACCESSION based fileNames
        "star_index":   {
                        "hg19": {
                                "female":   "hg19_female_v19_ERCC_sRNA_starIndex.tgz",
                                "male":     "hg19_male_v19_ERCC_sRNA_starIndex.tgz",
                                },
                        #"mm10": {
                        #        "female":   "mm10_male_sRNA_starIndex.tgz",
                        #        "male":     "mm10_male_sRNA_starIndex.tgz",
                        #        },
                        },
        "annotations":  {   # TODO: gender??  annotation??
                        "hg19": {
                                "v19":      "gencodeV19-tRNAs-ERCC.gtf.gz",
                                },
                        #"mm10": {
                        #        "M4":       "gencode.vM4-tRNAs-ERCC.gtf.gz",
                        #        },
                        },
        "chrom_sizes":  {
                        "hg19": {
                                "female":   "female.hg19.chrom.sizes",
                                "male":     "male.hg19.chrom.sizes",
                                },
                        "mm10": {
                                "female":   "male.mm10.chrom.sizes",
                                "male":     "male.mm10.chrom.sizes",
                                },
                        },
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

        # Now add pipline specific variables and tests
        
        # Could be multiple annotations supported per genome
        psv['annotation'] = args.annotation
        if psv['genome'] != self.GENOME_DEFAULT and psv['annotation'] == self.ANNO_DEFAULT:
            psv['annotation'] = self.ANNO_DEFAULTS[psv['genome']]
        if psv['annotation'] not in self.ANNO_ALLOWED[psv['genome']]:
            print psv['genome']+" has no "+psv['annotation']+" annotation."
            sys.exit(1)
        
        # Paired ends?
        if psv['paired_end']:
            print "Small-RNA is always expected to be single-end but mapping says otherwise."
            #print json.dumps(psv,indent=4,sort_keys=True)
            sys.exit(1)

        # Some specific settings
        psv['nthreads']   = 8

        # If annotation is not default, then add it to title
        if psv['annotation'] != self.ANNO_DEFAULTS[psv['genome']]:
            psv['title'] += ', ' + psv['annotation']
            psv['name']  += '_' + psv['annotation']

        # Must override results location because of annotation
        genome = psv['genome']
        if self.no_refs: # (no_refs is only True when templating)
            genome = None # If templating with no refs then this will hide genome and annotation
        psv['resultsLoc'] = self.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,psv['exp_type'], \
                                                                                        psv['genome'],psv['annotation'])
        psv['resultsFolder'] = psv['resultsLoc']
        if not self.template:
            psv['resultsFolder'] += psv['experiment'] + '/'
        self.update_rep_result_folders(psv)

        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        star_path = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][self.psv['gender']]
        star_fid = self.find_file(star_path,self.REF_PROJECT_DEFAULT)
        if star_fid == None:
            sys.exit("ERROR: Unable to locate STAR index file '" + star_path + "'")
        else:
            priors['star_index'] = star_fid

        anno_path = self.psv['refLoc']+self.REFERENCE_FILES['annotations'][self.psv['genome']][self.psv['annotation']]
        anno_fid = self.find_file(anno_path,self.REF_PROJECT_DEFAULT)
        if anno_fid == None:
            sys.exit("ERROR: Unable to locate Annotation file '" + anno_path + "'")
        else:
            priors['annotations'] = anno_fid

        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chrom_sizes_fid = self.find_file(chrom_sizes,self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        return priors

    #######################

if __name__ == '__main__':
    '''Run from the command line.'''
    srnaLaunch = SrnaLaunch()
    srnaLaunch.run()

