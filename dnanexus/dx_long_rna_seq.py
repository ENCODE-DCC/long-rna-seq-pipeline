import argparse
import os
import sys
import subprocess

import dxpy


ENCODE_DNA_ME_PROJECT_NAME = 'dna-me-pipeline'
''' This DNA Nexus project holds all the created applets and folders'''

ENCODE_REFERENCES_PROJECT = 'Encode Reference Files'
''' This DNA Nexus project holds Reference files (unaccessioned) used across several pipelines'''

ENCODE_SNAPSHOT_PROJECT = 'ENCODE-SDSC-snapshot-20140505'
''' This DNA Nexus project holds ENCFF files; should be replaced by a more permanent store '''

GENOME_REFERENCES = {
# Note this should be referred to by: biosample.donor.organism.name for any dataset
    'mouse':  'mm10.fa.gz',
    'human':  'hg19.fa.gz',
    'test':   'chr21.fa.gz'
}

REFERENCE_FILES = {}
APPLETS = {}

# TODO - load from pipeline object or .json text mockups
ANALYSIS_STEPS = [
    'index',
    'trim',
    'map',
    'extract'
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
                    default='human',
                    required=False)

    ap.add_argument('-t', '--test',
                    help='Use test input folder',
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

def get_project(project_name):
    project = dxpy.find_projects(name=project_name, name_mode='glob', return_handler=False)

    project = [p for p in project]
    if len(project) < 1:
        project = dxpy.DXProject(dxpy.api.project_new({'name': project_name, 'summary': 'ChIP-Seq Pipeline'})['id'])
    elif len(project) > 1:
        print 'Found more than 1 project matching ' + project_name + '.'
        print 'Please provide a unique project!'
        sys.exit(1)
    else:
        project = project[0]

    return project

def populate_workflow(wf, replicates, experiment, paired, gender, organism, applets_project_id):
    '''This function will populate the workflow for the methyl-seq Pipeline.'''

    genome = find_reference_file_by_name(GENOME_REFERENCES[organism], ENCODE_REFERENCES_PROJECT)
    # TODO somethink like loop over analysis_steps in pipeline objects
    ### INDEX
    index_input = {
        'genome': genome
    }
    stage_id = wf.add_stage(find_applet_by_name('index', applets_project_id), stage_input=index_input, folder=experiment)
    index_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'meIndex'
    })
    ### TRIM
    if not paired:
        trim_input = {
            'reads': replicates
        }

        stage_id = wf.add_stage(find_applet_by_name('trim-se', applets_project_id), stage_input=trim_input, folder=experiment)
        trim_output = dxpy.dxlink({
            'stage': stage_id,
            'outputField': 'trimmed_reads'
        })

        ### MAP
        map_input = {
            'genome': genome,
            'trimmed_reads': trim_output,
            'meIndex': index_output
        }
        stage_id = wf.add_stage(find_applet_by_name('map-se', applets_project_id), stage_input=map_input, folder=experiment)
        map_output = dxpy.dxlink({
            'stage': stage_id,
            'outputField': 'mapped_files'
        })
    else:
        if len(replicates) != 2:
            print "Must have exactly 2 replicats for paired-end pipeline"
            exit(1)

        trim_input = {
            'pair1_reads': replicates[0],
            'pair2_reads': replicates[1]
        }
        stage_id = wf.add_stage(find_applet_by_name('trim-pe', applets_project_id), stage_input=trim_input, folder=experiment)
        trim1_output = dxpy.dxlink({
            'stage': stage_id,
            'outputField': 'trimmed_reads1'
        })
        trim2_output = dxpy.dxlink({
            'stage': stage_id,
            'outputField': 'trimmed_reads2'
        })

        ### MAP
        map_input = {
            'genome': genome,
            'pair_1': trim1_output,
            'pair_2': trim2_output,
            'meIndex': index_output
        }
        stage_id = wf.add_stage(find_applet_by_name('map-pe', applets_project_id), stage_input=map_input, folder=experiment)
        map_output = dxpy.dxlink({
            'stage': stage_id,
            'outputField': 'mapped_files'
        })

    ### EXTRACT
    extract_input = {
        'genome': genome,
        'mapped_files': map_output
    }
    stage_id = wf.add_stage(find_applet_by_name('extract', applets_project_id), stage_input=extract_input, folder=experiment)


def copy_files(fids, project_id, folder):
    new_fids = []
    for file_dict in fids:
        (pid, fid) = file_dict.values()[0].values()
        f = dxpy.DXFile(dxid=fid,project=pid)
        fn = f.describe()['name']
        found_file = dxpy.find_one_data_object(classname='file', project=project_id, folder=folder, zero_ok=True, name=fn)

        if found_file is None:
            new_fids += [dxpy.dxlink(f.clone(project_id, folder))]
        else:
            new_fids += [dxpy.dxlink(found_file)]

    return new_fids

def project_has_folder(project, folder):
    folders = project.list_folder()['folders']

    return folder in folders

def resolve_applets_project():
    try:
        project = dxpy.find_one_project(name=ENCODE_DNA_ME_PROJECT_NAME, name_mode='exact', return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named {0}.'.format(ENCODE_DNA_ME_PROJECT_NAME)
        sys.exit(1)

    return dxpy.DXProject(project['id'])

def main():
    args = get_args()

    project = resolve_applets_project()
    #project = get_project(args.project_name)
    #print project.keys()
    print 'Project: ' + project.describe()['name']
    #print project.keys()
    print 'Experiment to analyze: ' + args.experiment
    project_folder = project_has_folder(project, '/'+args.experiment)
    if not project_folder:
        project_folder = project.new_folder('/'+args.experiment)

    #TODO get all replicate ids from encoded DB from ENCSR (args.experiment)
    #TODO error out if ENCSR not found, status not complete etc.
    if args.test:
        source_name = ENCODE_DNA_ME_PROJECT_NAME
        source_id = project.get_id()
    else:
        source_name = ENCODE_SNAPSHOT_PROJECT
        source_prj = dxpy.find_one_project(name=source_name, name_mode='exact', return_handler=False, level='VIEW')
        source_id = source_prj['id']


    if (len(args.replicates) < 1):
        sys.exit('Need to have at least 1 replicate file.')

    replicates = []
    for rep in args.replicates:
        dx_rep = dxpy.find_data_objects(classname='file', name=rep,
                                        name_mode='glob', project=source_id,
                                        return_handler=False)
        replicates.extend([ dxpy.dxlink(r) for r in dx_rep ])

    if not args.test:
        replicates = copy_files(replicates, project.id, "/"+args.experiment)

    if not replicates:
        print "No replicates found in project: " + project.name
        print "Looking for " + ", ".join(args.replicates)
        sys.exit(1)


    paired = args.paired
    gender = args.gender
    organism = 'human'
    #TODO determine paired or gender from ENCSR metadata
    # Now create a new workflow ()
    spec_name = args.experiment+'-'+'-'.join([ r.split('.')[0] for r in args.replicates])
    wf = dxpy.new_dxworkflow(title='dx_dna_me_'+spec_name,
                             name='ENCODE Bismark DNA-ME pipeline: '+spec_name,
                             description='The ENCODE Bismark pipeline for WGBS shotgun methylation analysis for experiment' + args.experiment,
                             folder='/'+args.experiment,
                             project=project.get_id())

    populate_workflow(wf, replicates, args.experiment, paired, gender, organism, project.id)
    #TODO - run the workflow automatically
    #TODO - export template workflows

if __name__ == '__main__':
    main()

