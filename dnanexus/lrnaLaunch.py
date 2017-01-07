#!/usr/bin/env python
# lrnaLaunch.py 2.1.1

import sys
import json

# import dxpy
from launch import Launch
# from template import Launch # (does not use dxencode at all)


class LrnaLaunch(Launch):
    '''Descendent from Launch class with 'long-rna-seq' methods'''

    PIPELINE_NAME = "long-rna-seq"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline analysis for one or more replicates. "
    ''' This pipline will compare only exactly two replicates replicates.'''

    GENOMES_SUPPORTED = ['hg19', 'GRCh38', 'mm10']
    ANNO_DEFAULTS = {'hg19': 'v19', 'GRCh38': 'v24', 'mm10': 'M4'}
    ANNO_ALLOWED = {'hg19':   [ANNO_DEFAULTS['hg19']],
                    'GRCh38': [ANNO_DEFAULTS['GRCh38']],
                    'mm10':   [ANNO_DEFAULTS['mm10'], 'M2', 'M3']}
    ANNO_DEFAULT = ANNO_DEFAULTS[Launch.GENOME_DEFAULT]
    ''' Multiple annotations might be supported for each genome.'''

    PIPELINE_BRANCH_ORDER = ["REP", "COMBINED_REPS"]
    '''Pipeline has standard replicate level processing and then combined replicate processing.'''

    PIPELINE_BRANCHES = {
        # '''Each branch must define the 'steps' and their (artificially) linear order.'''
        "REP": {
                "ORDER": {
                    "se":  ["align-tophat-se", "b2bw-se-top", "align-star-se", "b2bw-se-star",
                            "quant-rsem-alt"],
                    "pe":  ["align-tophat-pe", "b2bw-pe-top", "align-star-pe", "b2bw-pe-star", "quant-rsem"]
                },
                "STEPS": {
                            "align-tophat-se": {
                                        "app":     "align-tophat-se",
                                        "params":  {"library_id":   "library_id"},  # "nthreads"
                                        "inputs":  {"reads1":       "reads",
                                                    "tophat_index": "tophat_index"},
                                        "results": {"tophat_bam":   "tophat_bam"}
                            },
                            "align-tophat-pe": {
                                        "app":     "align-tophat-pe",
                                        "params":  {"library_id":   "library_id"},  # "nthreads"
                                        "inputs":  {"reads1":       "reads1",
                                                    "reads2":       "reads2",
                                                    "tophat_index": "tophat_index"},
                                        "results": {"tophat_bam":   "tophat_bam"}
                            },
                            "b2bw-se-top":  {
                                        "app":     "bam-to-bigwig-se-tophat",
                                        "params":  {"stranded":       "stranded"},
                                        "inputs":  {"tophat_bam":     "bam_file",
                                                    "chrom_sizes":    "chrom_sizes"},
                                        "results": {"tophat_all_bw":  "all_bw",
                                                    "tophat_uniq_bw": "uniq_bw"}
                            },
                            "b2bw-pe-top":  {
                                        "app":     "bam-to-bigwig-tophat",
                                        "params":  {"stranded":             "stranded"},
                                        "inputs":  {"tophat_bam":           "bam_file",
                                                    "chrom_sizes":          "chrom_sizes"},
                                        "results": {"tophat_minus_all_bw":  "minus_all_bw",
                                                    "tophat_minus_uniq_bw": "minus_uniq_bw",
                                                    "tophat_plus_all_bw":   "plus_all_bw",
                                                    "tophat_plus_uniq_bw":  "plus_uniq_bw"}
                            },
                            "align-star-se":   {
                                        "app":     "align-star-se",
                                        "params":  {"library_id":      "library_id"},  # "nthreads"
                                        "inputs":  {"reads1":          "reads",
                                                    "star_index":      "star_index"},
                                        "results": {"star_genome_bam": "star_genome_bam",
                                                    "star_anno_bam":   "star_anno_bam"}
                            },
                            "align-star-pe":   {
                                        "app":     "align-star-pe",
                                        "params":  {"library_id":      "library_id"},  # "nthreads"
                                        "inputs":  {"reads1":          "reads1",
                                                    "reads2":          "reads2",
                                                    "star_index":      "star_index"},
                                        "results": {"star_genome_bam": "star_genome_bam",
                                                    "star_anno_bam":   "star_anno_bam"}
                            },
                            "b2bw-se-star": {
                                        "app":     "bam-to-bigwig-se",
                                        "params":  {"stranded":       "stranded"},
                                        "inputs":  {"star_genome_bam": "bam_file",
                                                    "chrom_sizes":     "chrom_sizes"},
                                        "results": {"star_all_bw":     "all_bw",
                                                    "star_uniq_bw":    "uniq_bw"}
                            },
                            "b2bw-pe-star": {
                                        "app":     "bam-to-bigwig",
                                        "params":  {"stranded":           "stranded"},
                                        "inputs":  {"star_genome_bam":    "bam_file",
                                                    "chrom_sizes":        "chrom_sizes"},
                                        "results": {"star_minus_all_bw":  "minus_all_bw",
                                                    "star_minus_uniq_bw": "minus_uniq_bw",
                                                    "star_plus_all_bw":   "plus_all_bw",
                                                    "star_plus_uniq_bw":  "plus_uniq_bw"}
                            },
                            "quant-rsem":     {
                                        "app":     "quant-rsem",
                                        "params":  {"paired_end":       "paired_end",
                                                    "read_strand":       "read_strand"},
                                        "inputs":  {"star_anno_bam":     "star_anno_bam",
                                                    "rsem_index":        "rsem_index"},
                                        "results": {"rsem_iso_results":  "rsem_iso_results",
                                                    "rsem_gene_results": "rsem_gene_results"}
                            },
                            "quant-rsem-alt":     {
                                        "app":     "quant-rsem-alt",
                                        "params":  {"paired_end":       "paired_end",
                                                    "read_strand":       "read_strand"},
                                        "inputs":  {"star_anno_bam":     "star_anno_bam",
                                                    "rsem_index":        "rsem_index"},
                                        "results": {"rsem_iso_results":  "rsem_iso_results",
                                                    "rsem_gene_results": "rsem_gene_results"}
                            }
                }
        },
        "COMBINED_REPS": {
                "ORDER": {
                    "se": ["mad-qc-alt"],
                    "pe": ["mad-qc"]
                },
                "STEPS": {
                            "mad-qc": {
                                        "app":     "mad-qc",
                                        "params":  {},
                                        "inputs":  {"quants_a": "quants_a",
                                                    "quants_b": "quants_b"},
                                        "results": {"mad_plot": "mad_plot"}
                            },
                            "mad-qc-alt": {
                                        "app":     "mad-qc-alt",
                                        "params":  {},
                                        "inputs":  {"quants_a": "quants_a",
                                                    "quants_b": "quants_b"},
                                        "results": {"mad_plot": "mad_plot"}
                            }
                }
        }
    }

    PRUNE_STEPS = ["align-tophat-se", "align-tophat-pe","b2bw-se-top","b2bw-pe-top"]
    '''If --no-tophat is requested, these steps are pruned from the pipeline before launching.'''

    FILE_GLOBS = {
        # For looking up previous result files, use wild-cards
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
        "rsem_gene_results":    "/*_rsem.genes.results",
        "quants_a":             "/*_rsem.genes.results",
        "quants_b":             "/*_rsem.genes.results",
        "mad_plot":             "/*_mad_plot.png",
        }

    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should remove annotation if only one per genome
        # TODO: should use ACCESSION based fileNames
        "tophat_index":  {
                        "GRCh38": {
                                "female":   {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_tophatIndex.tgz"},
                                "male":     {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_tophatIndex.tgz"}
                                },
                        "hg19": {
                                "female":   {"v19": "hg19_female_v19_ERCC_tophatIndex.tgz"},
                                "male":     {"v19": "hg19_male_v19_ERCC_tophatIndex.tgz"}
                                },
                        "mm10": {
                                "female":   {
                                            "M2":  "mm10_male_M2_ERCC_tophatIndex.tgz",
                                            "M3":  "mm10_male_M3_ERCC_tophatIndex.tgz",
                                            "M4":  "mm10_XY_M4_ERCC_phiX_tophatIndex.tgz"
                                            },
                                "male":     {
                                            "M2":  "mm10_male_M2_ERCC_tophatIndex.tgz",
                                            "M3":  "mm10_male_M3_ERCC_tophatIndex.tgz",
                                            "M4":  "mm10_XY_M4_ERCC_phiX_tophatIndex.tgz"
                                            }
                                }
                        },
        "star_index":    {
                        "GRCh38": {
                                "female":   {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_starIndex.tgz"},
                                "male":     {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_starIndex.tgz"}
                                },
                        "hg19": {
                                "female":   {"v19": "hg19_female_v19_ERCC_starIndex.tgz"},
                                "male":     {"v19": "hg19_male_v19_ERCC_starIndex.tgz"}
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
        "rsem_index":    {
                        "GRCh38":   {"v24": "GRCh38_v24pri_tRNAs_ERCC_phiX_rsemIndex.tgz"},
                        "hg19":     {"v19": "hg19_male_v19_ERCC_rsemIndex.tgz"},
                        "mm10":     {
                                    "M2":  "mm10_male_M2_ERCC_rsemIndex.tgz",
                                    "M3":  "mm10_male_M3_ERCC_rsemIndex.tgz",
                                    "M4":  "mm10_XY_M4_ERCC_phiX_rsemIndex.tgz"
                                    }
                        },
        "chrom_sizes":   {
                        "GRCh38":   {"female":   "GRCh38_EBV.chrom.sizes",
                                     "male":     "GRCh38_EBV.chrom.sizes"},
                        "hg19":     {"female":   "female.hg19.chrom.sizes",
                                     "male":     "male.hg19.chrom.sizes"},
                        "mm10":     {"female":   "mm10_no_alt.chrom.sizes",
                                     "male":     "mm10_no_alt.chrom.sizes"}
                        }
        }

    def __init__(self):
        Launch.__init__(self)

    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self, parse=False)

        ap.add_argument('-a', '--annotation',
                        help="Label of annotation (default: '" + self.ANNO_DEFAULT + "')",
                        choices=[self.ANNO_DEFAULT, 'M2', 'M3', 'M4'],
                        default=self.ANNO_DEFAULT,
                        required=False)

        ap.add_argument('--no_tophat',  # This has become a noop retained for consistency with hg19
                        help='Do not include TopHat steps in pipeline (default: exclude TopHat steps).',
                        action='store_true',
                        required=False)

        ap.add_argument('--tophat_also',
                        help='Do not include TopHat steps in pipeline (default: exclude TopHat steps).',
                        action='store_true',
                        required=False)

        return ap.parse_args()

    def pipeline_specific_vars(self, args, verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self, args)

        # Could be multiple annotations supported per genome
        psv['annotation'] = args.annotation
        if psv['genome'] != self.GENOME_DEFAULT and psv['annotation'] == self.ANNO_DEFAULT:
            psv['annotation'] = self.ANNO_DEFAULTS[psv['genome']]
        if psv['annotation'] not in self.ANNO_ALLOWED[psv['genome']]:
            print psv['genome']+" has no "+psv['annotation']+" annotation."
            sys.exit(1)

        # Some specific settings
        psv['nthreads'] = 8
        psv['rnd_seed'] = 12345

        # If paired-end then read_strand might vary TruSeq or ScriptSeq, but only for quant-rsem
        psv["read_strand"] = "unstranded"  # SE experiments are all unstranded
        if psv["paired_end"]:
            psv["read_strand"] = "reverse"  # Usual ENCODE LRNA experiments are rd1-/rd2+ (AKA reverse)
            if not psv["stranded"]:
                psv["read_strand"] = "unstranded"  # "ScriptSeq" experiments are rd1+/rd2- (AKA forward)
                print "Detected unstranded library"
            elif psv.get('ScriptSeq', False):  # file.replicate.library.document contains "/documents/F17c31e10-1542-42c6-8b4c-3afff95564cf%2F"
                psv["read_strand"] = "ScriptSeq"  # "ScriptSeq" experiments are rd1+/rd2- (AKA forward)
                print "Detected ScriptSeq"
        # print "Detected special cases"

        # If annotation is not default, then add it to title
        if psv['annotation'] != self.ANNO_DEFAULTS[psv['genome']]:
            psv['title'] += ', ' + psv['annotation']
            psv['name'] += '_' + psv['annotation']

        self.no_tophat = True
        if args.tophat_also:
            self.no_tophat = False
            self.PRUNE_STEPS = []  # This blocks pruning... keeping tophat

        # Must override results location because of annotation
        psv['resultsLoc'] = self.umbrella_folder(args.folder, self.FOLDER_DEFAULT, self.proj_name,
                                                 psv['exp_type'], psv['genome'], psv['annotation'])
        psv['resultsFolder'] = psv['resultsLoc']
        if not self.template:
            psv['resultsFolder'] += psv['experiment'] + '/'
        self.update_rep_result_folders(psv)

        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv, indent=4)
        return psv

    def find_ref_files(self, priors):
        '''Locates all reference files based upon gender, organism and annotation.'''
        top_path = self.psv['refLoc']+self.REFERENCE_FILES['tophat_index'][self.psv['genome']][
                                                    self.psv['gender']][self.psv['annotation']]
        top_fid = self.find_file(top_path, self.REF_PROJECT_DEFAULT)
        if top_fid is None:
            sys.exit("ERROR: Unable to locate TopHat index file '" + top_path + "'")
        else:
            priors['tophat_index'] = top_fid

        star_path = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][
                                                    self.psv['gender']][self.psv['annotation']]
        star_fid = self.find_file(star_path, self.REF_PROJECT_DEFAULT)
        if star_fid is None:
            sys.exit("ERROR: Unable to locate STAR index file '" + star_path + "'")
        else:
            priors['star_index'] = star_fid

        rsem_path = self.psv['refLoc']+self.REFERENCE_FILES['rsem_index'][self.psv['genome']][
                                                                        self.psv['annotation']]
        rsem_fid = self.find_file(rsem_path, self.REF_PROJECT_DEFAULT)
        if rsem_fid is None:
            sys.exit("ERROR: Unable to locate RSEM index file '" + rsem_path + "'")
        else:
            priors['rsem_index'] = rsem_fid

        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][
                                                                                self.psv['gender']]
        chrom_sizes_fid = self.find_file(chrom_sizes, self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid is None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        return priors

    # ######################


if __name__ == '__main__':
    '''Run from the command line.'''
    lrnaLaunch = LrnaLaunch()
    lrnaLaunch.run()

