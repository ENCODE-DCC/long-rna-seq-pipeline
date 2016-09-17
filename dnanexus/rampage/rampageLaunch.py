#!/usr/bin/env python
# rampageLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
from launch import Launch
#from template import Launch # (does not use dxencode at all)

class RampageLaunch(Launch):
    '''Descendent from Launch class with 'rampage' methods'''

    PIPELINE_NAME = "rampage"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''
    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    GENOMES_SUPPORTED = ['hg19', 'GRCh38', 'mm10']
    ANNO_DEFAULTS = {'hg19': 'v19', 'GRCh38': 'v24', 'mm10': 'M4' }
    ANNO_ALLOWED = { 'hg19': [ ANNO_DEFAULTS['hg19'] ],
                     'GRCh38': [ ANNO_DEFAULTS['GRCh38'] ],
                     'mm10': [ ANNO_DEFAULTS['mm10'], 'M2', 'M3' ] }
    ANNO_DEFAULT = ANNO_DEFAULTS[Launch.GENOME_DEFAULT]
    ''' Multiple annotations might be supported for each genome.'''

    RESULT_FOLDER_DEFAULT = '/rampage/'
    ''' This the default location to place results folders for each experiment.'''

    # NOTE '/lrna/' or '/run/' is safer and more efficient than '/', but also more brittle
    #CONTROL_ROOT_FOLDER = '/'
    CONTROL_ROOT_FOLDER = '/long-RNA-seq/'
    ''' Rampage requires a control file which may be discoverable.'''
    CONTROL_FILE_GLOB = '*_star_genome.bam'

    PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
    '''This pipeline has the standard replicate level processing and then combined replicate processing.'''

    PIPELINE_BRANCHES = {
    #'''Each branch must define the 'steps' and their (artificially) linear order.'''
         "REP": {
                "ORDER": [ "rampage-align-pe", "rampage-signals", "rampage-peaks" ],
                "STEPS": {
                            "rampage-align-pe": {
                                "app":     "rampage-align-pe",
                                # lrna_221 no assay "params":  { "library_id": "library_id", "assay_type": "assay_type", "nthreads": "nthreads" },
                                "params":  { "library_id": "library_id", "nthreads": "nthreads" },
                                "inputs":  { "reads1":     "reads1",
                                             "reads2":     "reads2",
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
                                             "unique_minus_bw": "unique_minus_bw" }
                            },
                            "rampage-peaks": {
                                "app":     "rampage-peaks",
                                # lrna_221 no assay "params":  { "assay_type": "assay_type", "nthreads": "nthreads" },
                                "params":  { "nthreads": "nthreads" },
                                "inputs":  { "control_bam": "control_bam", 
                                             "gene_annotation": "gene_annotation", 
                                             "chrom_sizes": "chrom_sizes",
                                             "rampage_marked_bam": "rampage_marked_bam" },
                                "results": { "rampage_peaks_bed":  "rampage_peaks_bed",
                                             "rampage_peaks_bb":   "rampage_peaks_bb",
                                             "rampage_peaks_gff":  "rampage_peaks_gff",
                                             "rampage_peak_quants":"rampage_peak_quants" }
                            }
                }
         },
         "COMBINED_REPS": {
                "ORDER": [ "rampage-idr", "rampage-mad-qc" ],
                "STEPS": {
                            "rampage-idr": {
                                "app":     "rampage-idr",
                                "params":  {},
                                "inputs":  { "peaks_a": "peaks_a", 
                                             "peaks_b": "peaks_b", 
                                             "chrom_sizes": "chrom_sizes" },
                                "results": { "rampage_idr_png": "rampage_idr_png",
                                             "rampage_idr_bb":  "rampage_idr_bb",
                                             "rampage_idr_bed": "rampage_idr_bed" }
                            },
                            "rampage-mad-qc": {
                                        "app":     "rampage-mad-qc",
                                        "params":  {},
                                        "inputs":  { "quants_a": "quants_a", 
                                                     "quants_b": "quants_b" },
                                        "results": { "mad_plot": "mad_plot" }
                            },
                }
         }
    }

    FILE_GLOBS = {
        "all_plus_bw":        "/*_5p_plusAll.bw",
        "rampage_marked_bam": "/*_star_marked.bam",
        "all_minus_bg":       "/*_5p_minusAll.bg",
        "unique_plus_bg":     "/*_5p_plusUniq.bg",
        "rampage_peaks_bed":  "/*_peaks.bed.gz",
        "rampage_peaks_bb":   "/*_peaks.bb",
        "rampage_peaks_gff":  "/*_peaks.gff.gz",
        "rampage_peak_quants":"/*_peaks_quant.tsv",
        "peaks_a":            "/*_peaks.bed.gz",
        "peaks_b":            "/*_peaks.bed.gz",
        "quants_a":           "/*_peaks_quant.tsv",
        "quants_b":           "/*_peaks_quant.tsv",
        "all_plus_bg":        "/*_5p_plusAll.bg",
        "rampage_star_log":   "/*_star_Log.final.out",
        "all_minus_bw":       "/*_5p_minusAll.bw",
        "unique_plus_bw":     "/*_5p_plusUniq.bw",
        "unique_minus_bg":    "/*_5p_minusUniq.bg",
        "unique_minus_bw":    "/*_5p_minusUniq.bw",
        "rampage_idr_bed":    "/*_idr.bed.gz",
        "rampage_idr_bb":     "/*_idr.bb",
        "rampage_idr_png":    "/*_idr.png",
        "mad_plot":           "/*_mad_plot.png",
    }

    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should remove annotation if only one per genome
        # TODO: should use ACCESSION based fileNames
        "star_index":   {
                        "GRCh38": {
                                "female":   {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_starIndex.tgz"},
                                "male":     {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_starIndex.tgz"}
                                },
                        "hg19": {
                                "female":   {"v19": "hg19_female_v19_ERCC_starIndex.tgz"         },
                                "male":     {"v19": "hg19_male_v19_ERCC_starIndex.tgz"           }
                                },
                        "mm10": {
                                "female":   {
                                            "M2":  "mm10_male_M2_ERCC_starIndex.tgz",
                                            "M3":  "mm10_male_M3_ERCC_starIndex.tgz",
                                            "M4":  "mm10_XY_M4_ERCC_phiX_starIndex.tgz"
                                            },
                                "male":     {
                                            "M2":  "mm10_male_M2_ERCC_starIndex.tgz",
                                            "M3":  "mm10_male_M3_ERCC_starIndex.tgz",
                                            "M4":  "mm10_XY_M4_ERCC_phiX_starIndex.tgz"
                                            }
                                }
                        },
        "gene_annotation":   {
                        "GRCh38":   {"v24": "gencode.v24.primary_assembly.annotation.gtf.gz"},
                        "hg19":     {"v19": "gencode.v19.annotation.gtf.gz"                 },
                        "mm10":     {
                                    "M2":  "gencode.vM2.annotation.gtf.gz",
                                    "M3":  "gencode.vM3.annotation.gtf.gz",
                                    "M4":  "gencode.vM4.annotation.gtf.gz"
                                    }
                         },
        "chrom_sizes":   {
                        "GRCh38":   {"female":   "GRCh38_EBV.chrom.sizes",
                                     "male":     "GRCh38_EBV.chrom.sizes"  },
                        "hg19":     {"female":   "female.hg19.chrom.sizes",
                                     "male":     "male.hg19.chrom.sizes"   },
                        "mm10":     {"female":   "mm10_no_alt.chrom.sizes",
                                     "male":     "mm10_no_alt.chrom.sizes"   }
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
        #args.pe = True # This is necessary to ensure templating does what it must.
        psv = Launch.pipeline_specific_vars(self,args)
        
        # Could be multiple annotations supported per genome
        psv['annotation'] = args.annotation
        if psv['genome'] != self.GENOME_DEFAULT and psv['annotation'] == self.ANNO_DEFAULT:
            psv['annotation'] = self.ANNO_DEFAULTS[psv['genome']]
        if psv['annotation'] not in self.ANNO_ALLOWED[psv['genome']]:
            print psv['genome']+" has no "+psv['annotation']+" annotation."
            sys.exit(1)

        # Some specific settings
        psv['assay_type'] = "rampage"
        if self.exp["assay_term_name"] == "CAGE":
            psv['assay_type'] = "cage"
        psv['nthreads']   = 8
        if not self.template:
            psv['control'] = args.control
        
        if psv['paired_end'] and psv['assay_type'] == "cage":
            print "CAGE is always expected to be single-end but mapping says otherwise."
            sys.exit(1)
        elif not psv['paired_end'] and psv['assay_type'] == "rampage":
            print "Rampage is always expected to be paired-end but mapping says otherwise."
            sys.exit(1)

        # run will either be for combined or single rep.
        if not self.combined_reps:
            run = psv['reps']['a']  # If not combined then run will be for the first (only) replicate
        else:
            run = psv
            
        # If annotation is not default, then add it to title
        if psv['annotation'] != self.ANNO_DEFAULTS[psv['genome']]:
            psv['title'] += ', ' + psv['annotation']
            psv['name']  += '_' + psv['annotation']

        if self.exp["assay_term_name"] == "CAGE":
            psv['name'] = psv['assay_type'] + psv['name'][4:]
            psv['title'] = "CAGE" + psv['title'][7:]

        # Must override results location because of annotation
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
        star_path = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][self.psv['gender']][self.psv['annotation']]
        star_fid = self.find_file(star_path,self.REF_PROJECT_DEFAULT)
        if star_fid == None:
            sys.exit("ERROR: Unable to locate STAR index file '" + star_path + "'")
        else:
            priors['star_index'] = star_fid

        anno_path = self.psv['refLoc']+self.REFERENCE_FILES['gene_annotation'][self.psv['genome']][self.psv['annotation']]
        anno_fid = self.find_file(anno_path,self.REF_PROJECT_DEFAULT)
        if anno_fid == None:
            sys.exit("ERROR: Unable to locate Gene Annotation file '" + anno_path + "'")
        else:
            priors['gene_annotation'] = anno_fid

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
    rampageLaunch = RampageLaunch()
    rampageLaunch.run()

