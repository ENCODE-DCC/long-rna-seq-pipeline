#!/usr/bin/env python
import argparse
import os
import sys
import subprocess

import dxpy


ENCODE_DNA_ME_PROJECT_NAME = 'long-rna-seq-pipeline'
''' This DNA Nexus project holds all the created applets and folders'''

ENCODE_REFERENCES_PROJECT = 'ENCODE Reference Files'
''' This DNA Nexus project holds Reference files (unaccessioned) used across several pipelines'''

ENCODE_SNAPSHOT_PROJECT = 'ENCODE-SDSC-snapshot-20140505'
''' This DNA Nexus project holds ENCFF files; should be replaced by a more permanent store '''

ENCODE_PUBLIC_PROJECT = 'ENCODE Universal Processing Pipelines'
PUBLIC_FOLDER = '/RNA-Seq (long)'

GENOME_REFERENCES = {
# Note this should be referred to by: biosample.donor.organism.name for any dataset
    'mouse':  {
        'm': {
            'genome': 'male.mm9.fa.gz',
            'gene_annotation': '',
            'trna_annotation': ''
        },
        'f': {
            'genome': 'female.mm9.fa.gz',
            'gene_annotation': '',
            'trna_annotation': ''
        }
    },
    'human':  {
        'm': {
            'genome': 'male.hg19.fa.gz',
            'gene_annotation': 'gencode.v19.annotation.gtf.gz',
            'trna_annotation': 'gencode.v19.tRNAs.gtf.gz'
        },
        'f': {
            'genome': 'female.hg19.fa.gz',
            'gene_annotation': 'gencode.v19.annotation.gtf.gz',
            'trna_annotation': 'gencode.v19.tRNAs.gtf.gz'
        }
    },
    'test':   'chr21.fa.gz'
}

REFERENCE_FILES = {}
APPLETS = {}

# TODO - load from pipeline object or .json text mockups
ANALYSIS_STEPS = [
    'merge-annotation',
    'prep-star',
    'prep-rsem',
    'prep-tophat',
    'align-star-se',
    'align-star-pe',
    'align-tophat-pe',
    'align-tophat-se',
    'quant-rsem',
    'bam-to-bigwig-stranded',
    'bam-to-bigwig-unstranded'
 ]

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Generate DNAnexus workflow for the ENCODE methylation pipeline.')

    ap.add_argument('-r', '--replicates',
                    help='Replicate fastq files.',
                    nargs='+',
                    required=False)

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    required=True)

    ap.add_argument('-g', '--gender',
                    help='Gender of sample',
                    choices=['m', 'f'],
                    default='m',
                    required=False)

    ap.add_argument('-p', '--paired',
                    help='Force use of paired-end pipeline',
                    action='store_true',
                    required=False)

    ap.add_argument('-o', '--organism',
                    help='Organism to map to',
                    required=True)

    ap.add_argument('-l', '--library',
                    help='ENCODE accession of biosample library (for BAM header)',
                    required=True)

    ap.add_argument('-t', '--test',
                    help='Use test input folder',
                    action='store_true',
                    required=False)

    ap.add_argument('-x', '--export',
                    help='Export generic Workflow (no inputs) to DNA Nexus project',
                    action='store_true',
                    required=False)

    return ap.parse_args()

def find_reference_file_by_name(reference_name, project_name):
    '''Looks up a reference file by name in the project that holds common tools. From Joe Dale's code.'''
    project = dxpy.find_one_project(name=project_name, name_mode='exact', return_handler=False)
    cached = '*'
    if (reference_name, project['id']) not in REFERENCE_FILES:
        found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                          project=project['id'],
                                          recurse=True,
                                          zero_ok=False, more_ok=False, return_handler=True)
        REFERENCE_FILES[(reference_name, project['id'])] = found
        cached = ''

    print cached + "Resolved %s to %s" % (reference_name, REFERENCE_FILES[(reference_name, project['id'])].get_id())
    return dxpy.dxlink(REFERENCE_FILES[(reference_name, project['id'])])


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '*'
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                          project=applets_project_id,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    print cached + "Resolved %s to %s" % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id())
    return APPLETS[(applet_name, applets_project_id)]


def populate_workflow(wf, replicates, experiment, inputs, applets_project_id, export):
    '''This function will populate the workflow for the methyl-seq Pipeline.'''

    gene_annotation = find_reference_file_by_name(GENOME_REFERENCES[inputs['organism']][inputs['gender']]['gene_annotation'], ENCODE_REFERENCES_PROJECT)
    trna_annotation = find_reference_file_by_name(GENOME_REFERENCES[inputs['organism']][inputs['gender']]['trna_annotation'], ENCODE_REFERENCES_PROJECT)
    spike_in = find_reference_file_by_name('ENCFF001RTP.fasta', ENCODE_DNA_ME_PROJECT_NAME) ### note not in reference project
    genome = find_reference_file_by_name(GENOME_REFERENCES[inputs['organism']][inputs['gender']]['genome'], ENCODE_REFERENCES_PROJECT)
    index_prefix = inputs['spec_name']

    #'merge-annotation',
    if export:
        merge_input = {}
    else:
        merge_input = {
            'gene_annotation': gene_annotation,
            'trna_annotation': trna_annotation,
            'spike_in': spike_in
        }
    stage_id = wf.add_stage(find_applet_by_name('merge-annotation', applets_project_id), stage_input=merge_input, folder=experiment)
    merge_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'combined_gtf'
    })
    #'prep-star', if export:
    prep_input = {
        'annotations': merge_output
    }
    if not export:
        prep_input['genome'] = genome
        prep_input['spike_in'] = spike_in
        prep_input['index_prefix'] = index_prefix
    stage_id = wf.add_stage(find_applet_by_name('prep-star', applets_project_id), stage_input=prep_input, folder=experiment)
    prep_star_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'star_index'
    })
    #'prep-rsem',
    stage_id = wf.add_stage(find_applet_by_name('prep-rsem', applets_project_id), stage_input=prep_input, folder=experiment)
    prep_rsem_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'rsem_index'
    })
    #'prep-tophat',
    stage_id = wf.add_stage(find_applet_by_name('prep-tophat', applets_project_id), stage_input=prep_input, folder=experiment)
    prep_tophat_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'tophat_index'
    })

    ## alignment steps
    align_input = {
    }
    align_input['library_id'] = inputs['library_id']
    if not inputs['paired']:
        star_step = 'align-star-se'
        th_step = 'align-tophat-se'
        align_input['reads'] = dxpy.dxlink(replicates[0])
    else:
        star_step = 'align-star-pe'
        th_step = 'align-tophat-pe'
        align_input['reads_1'] = dxpy.dxlink(replicates[0])
        align_input['reads_2'] = dxpy.dxlink(replicates[1])


    align_input['star_index'] = prep_star_output
    if export:
        del align_input['library_id']
        if align_input.get('reads'):
            del align_input['reads']
        elif align_input.get('reads_1'):
            del align_input['reads_1']
            del align_input['reads_2']
    stage_id = wf.add_stage(find_applet_by_name(star_step, applets_project_id), stage_input=align_input, folder=experiment)
    star_annotation_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'annotation_bam'
    })

    del align_input['star_index']
    align_input['tophat_index'] = prep_tophat_output
    stage_id = wf.add_stage(find_applet_by_name(th_step, applets_project_id), stage_input=align_input, folder=experiment)
    tophat_genome_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'genome_bam'
    })
    ## above file is needed for bam-to-bigwig (unimplemented)

    #'quant-rsem'
    quant_input = {
        'rsem_index': prep_rsem_output,
        'annotation_bam': star_annotation_output,
        'paired': inputs['paired'],
        'stranded': inputs['stranded']
    }
    if not export:
        quant_input['read_prefix'] = index_prefix
        quant_input['rnd_seed'] = inputs['rnd_seed']

    stage_id = wf.add_stage(find_applet_by_name('quant-rsem', applets_project_id), stage_input=quant_input, folder=experiment)

    ## below unimplemented
    #'bam-to-bigwig-stranded',
    #'bam-to-bigwig-unstranded'

def copy_files(fids, project_id, folder):
    new_fids = []
    for file_dict in fids:
        f = dxpy.DXFile(dxid=file_dict['id'], project=file_dict['project'])
        fn = f.describe()['name']

        # Check to see if file already exists.
        found_file = dxpy.find_one_data_object(classname='file', project=project_id, folder=folder, zero_ok=True, name=fn)
        if found_file is None:
            new_fids += [dxpy.dxlink(f.clone(project_id, folder))]
        else:
            new_fids += [dxpy.dxlink(found_file)]

    return new_fids

def project_has_folder(project, folder):
    folders = project.list_folder()['folders']

    return folder in folders

def resolve_project(project_name, level=None):
    try:
        project = dxpy.find_one_project(name=project_name, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named {0}.'.format(project_name)
        sys.exit(1)

    return dxpy.DXProject(project['id'])

def main():
    args = get_args()
    if len(args.replicates) < 1:
        sys.exit('Need to have at least 1 replicate file.')

    project = resolve_project(ENCODE_DNA_ME_PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    print 'Experiment to analyze: ' + args.experiment
    if not project_has_folder(project, '/'+args.experiment):
        project.new_folder('/'+args.experiment)

    #TODO get all replicate ids from encoded DB from ENCSR (args.experiment)
    #TODO error out if ENCSR not found, status not complete etc.
    if args.test:
        source_id = project.get_id()
    else:
        source_id = resolve_project(ENCODE_SNAPSHOT_PROJECT, level='VIEW').get_id()

    replicates = []
    for rep in args.replicates:
        dx_rep = dxpy.find_data_objects(classname='file', name=rep,
                                        name_mode='exact', project=source_id,
                                        return_handler=False)
        replicates.extend(dx_rep)

    if not args.test:
        replicates = copy_files(replicates, project.get_id(), "/"+args.experiment)

    if not replicates:
        print "No replicates found in project: " + project.name
        print "Looking for " + ", ".join(args.replicates)
        sys.exit(1)

    inputs = {
        'rnd_seed': 12345
    }
    inputs['paired'] = args.paired
    inputs['gender']= args.gender
    inputs['organism'] = args.organism
    inputs['library_id'] = args.library
    #TODO determine paired or gender from ENCSR metadata
    # Now create a new workflow ()
    inputs['spec_name'] = args.experiment+'-'+'-'.join([ r.split('.')[0] for r in args.replicates])
    title_root = 'dx_long_rna_seq_'
    name_root = 'ENCODE Long RNA Seq: '
    desc = 'The ENCODE RNA Seq pipeline for long RNAs'
    if args.paired:
        title_root = title_root + '_paired_end '
        name_root = name_root + '(paired-end) '
        inputs['stranded'] = True
    else:
        title_root = title_root + '_single_end '
        name_root = name_root + '(single-end) '
        inputs['stranded'] = False


    if args.export:
        project_id = dxpy.find_one_project(name=ENCODE_PUBLIC_PROJECT, name_mode='exact', return_handler=False)['id']
        wf = dxpy.new_dxworkflow(title=title_root,
                                 name=name_root,
                                 description=desc,
                                 folder=PUBLIC_FOLDER,
                                 project=project_id)
    else:
        project_id = project.get_id()
        wf = dxpy.new_dxworkflow(title=title_root+inputs['spec_name'],
                             name=name_root+inputs['spec_name'],
                             description=desc+' for experiment:' + args.experiment,
                             folder='/'+args.experiment,
                             project=project.get_id())

    populate_workflow(wf, replicates, args.experiment, inputs, project.id, args.export)
    #TODO - run the workflow automatically
    #TODO - export template workflows

if __name__ == '__main__':
    main()

