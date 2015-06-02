#!/usr/bin/env python
# srnaLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
#from dxencode import dxencode as dxencode
import dxencode as dxencode
from launch import Launch

class SrnaLaunch(Launch):
    '''Descendent from Launch class with 'small-rna-seq' methods'''

    PIPELINE_NAME = "small-rna-seq"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline analysis for one replicate. "
    ''' This pipline does not support combined replicates.'''
                    
    PIPELINE_BRANCH_ORDER = [ "REP" ]
    '''Currently there are no combined replicate processing steps.'''
    
    PIPELINE_BRANCHES = {
    #'''Each branch must define the 'steps' and their (artificially) linear order.'''
         "REP": {
                "ORDER": [ "concat-fastqs", "small-rna-align", "small-rna-signals" ]
                # "STEPS" can be discovered in this simple pipeline
         }
    }

    REFERENCE_FILES = {
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

    def __init__(self):
        Launch.__init__(self)
        
    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self,args)

        # Now add pipline specific variables and tests
        
        # Paired ends?
        if psv['paired_end']:
            print "Small-RNA is always expected to be single-end but mapping says otherwise."
            #print json.dumps(psv,indent=4,sort_keys=True)
            sys.exit(1)

        # Some specific settings
        psv['nthreads']   = 8

        # run will either be for combined or single rep.
        if self.combined_reps:
            print "Small-RNA-seq pipeline currently does not support combined-replicate processing."
            sys.exit(1)
           
        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        star_ix = self.psv['refLoc']+self.REFERENCE_FILES['star_index'][self.psv['genome']][self.psv['gender']]
        star_ix_fid = dxencode.find_file(star_ix,dxencode.REF_PROJECT_DEFAULT)
        if star_ix_fid == None:
            sys.exit("ERROR: Unable to locate STAR index file '" + star_ix + "'")
        else:
            priors['star_index'] = star_ix_fid

        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chrom_sizes_fid = dxencode.find_file(chrom_sizes,dxencode.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()

    #######################

if __name__ == '__main__':
    '''Run from the command line.'''
    srnaLaunch = SrnaLaunch()
    srnaLaunch.run()

