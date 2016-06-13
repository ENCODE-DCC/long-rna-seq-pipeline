#!/usr/bin/env python2.7
# parse_property.py  Reads property string and prses a requested value.
#                    Write request to stdout and verbose info to stderr.  This allows easy use in dx app scripts.

# imports needed for Settings class:
import os, sys, string, argparse, json
import dxpy

def env_get_current_project_id():
    ''' Returns the current project name for the command-line environment '''
    err, proj_name = commands.getstatusoutput('cat ~/.dnanexus_config/DX_PROJECT_CONTEXT_NAME')
    if err != 0:
        return None
    proj = dxencode.get_project(proj_name)
    return proj.get_id()
    
    return proj_name

def get_dxfile(filePath,project=None):
    '''Returns dxfile object.'''
    dxfile = None
    #if filePath.find("$dnanexus_link") != -1:
    #    filePath = filePath.split(' ')[1]
    #    filePath = filePath.replace("'","").replace('"','').replace("}","").replace("{","")
    try:
        dxlink = json.loads(filePath.strip("'"))
    except:
        dxlink = None
        
    if project != None:
        
        try:
            if dxlink != None:
                dxfile = dxpy.get_handler(dxlink,project=project)
            else:
                dxfile = dxpy.get_handler(filePath,project=project)
        except:
            try:
                dxlink = dxpy.dxlink(filePath,project=project)
                dxfile = dxpy.get_handler(dxlink)
            except:
                try:
                    proj_id = env_get_current_project_id()
                    dxfile = dxpy.DXFile(filePath,project=proj_id)
                except:
                    sys.stderr.write('ERROR: unable to find file "' + filePath + '": \n')
                    sys.exit(0)  # Do not error on tool run in dx script 
    
    else:
    
        try:
            if dxlink != None:
                dxfile = dxpy.get_handler(dxlink)
            else:
                dxfile = dxpy.get_handler(filePath)
        except:
            try:
                dxlink = dxpy.dxlink(filePath)
                dxfile = dxpy.get_handler(dxlink)
            except:
                try:
                    proj_id = env_get_current_project_id()
                    dxfile = dxpy.DXFile(filePath,project=proj_id)
                except:
                    sys.stderr.write('ERROR: unable to find file "' + filePath + '": \n')
                    sys.exit(0)  # Do not error on tool run in dx script 

    if dxfile == None:
        sys.stderr.write('ERROR: unable to find file "' + filePath + '": \n')
        sys.exit(0)  # Do not error on tool run in dx script 
    
    return dxfile


def file_get_property(filePath,key,subkey,return_json=False,project=None,verbose=False):
    '''Returns dx file's property matching 'key'.'''

    dxfile = get_dxfile(filePath,project=project)
    
    props = dxfile.get_properties()
    if not props:
        sys.stderr.write('ERROR: unable to find properties for file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script 
    
    if key not in props:
        sys.stderr.write('ERROR: unable to find "'+key+'" in properties for file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script
    props = props[key]
         
    if return_json or subkey != None:
        try:
            props = json.loads(props)
        except:
            try:
                props = json.loads("{"+props+"}")
            except:
                sys.stderr.write('Failure parsing "'+props+'" as json.\n') 
                sys.exit(0)  # Do not error on tool run in dx script

    if subkey != None:
        if subkey not in props:
            sys.stderr.write('ERROR: unable to find "'+subkey+'" in properties for file "' + filePath + '": \n') 
            sys.exit(0)  # Do not error on tool run in dx script
        props = props[subkey]
        
    if verbose:
        sys.stderr.write(props + '\n')
    
    return props

def file_describe(filePath,key=None,project=None,verbose=False):
    '''Returns dx file's description property matching 'key'.'''

    dxfile = get_dxfile(filePath,project=project)    
    
    desciption = dxfile.describe()
    if not desciption:
        sys.stderr.write('ERROR: unable to find description of file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script 
    
    if key == None:
        if verbose:
            sys.stderr.write(json.dumps(desciption) + '\n')
        return desciption
    
    if key not in desciption:
        sys.stderr.write('ERROR: unable to find "'+key+'" in description of file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script
    value = desciption[key]
         
    if verbose:
        sys.stderr.write(value + '\n')
    
    return value

def file_details(filePath,key=None,project=None,verbose=False):
    '''Returns dx file's description property matching 'key'.'''

    dxfile = get_dxfile(filePath,project=project)    
    
    details = dxfile.describe(incl_details=True).get('details')
    if not details:
        sys.stderr.write('ERROR: unable to find details of file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script 
    
    if key == None:
        if verbose:
            sys.stderr.write(json.dumps(details) + '\n')
        return details
    
    if key not in details:
        sys.stderr.write('ERROR: unable to find "'+key+'" in details of file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script
    value = details[key]
         
    if verbose:
        sys.stderr.write(value + '\n')
    
    return value

def folder_create_root(folder,verbose=False):
    '''Returns a standard file name root when in folder {exp_acc}/repN_N.'''

    if verbose:
        sys.stderr.write("Trying to create root_name from: "+ folder + '\n')
    
    folders = folder.split('/')
    exp = ''
    rep = ''
    for one_folder in reversed(folders):
        if one_folder.count(' ') > 0:
            continue
        if rep == '' and one_folder.startswith('rep'):
            rep = one_folder
        elif exp == '' and one_folder.startswith('ENCSR'):
            exp = one_folder
            break
    
    root = ''
    if len(exp) > 0:
        root = exp
        if len(rep) > 0:
            root += '_' + rep

    if verbose:
        sys.stderr.write(root + '\n')
    
    return root

def file_create_root(filePath,project=None,verbose=False,quiet=False):
    '''Returns a standard file name root when in folder {exp_acc}/repN_N.'''

    folder = file_describe(filePath,'folder',project=project,verbose=verbose)
    root = folder_create_root(folder,verbose=verbose)

    if root == '' and not quiet:  # Note, this is not an error when file is from a different dx context!
        sys.stderr.write("Found nothing for root: folder["+folder+"] path ["+filePath+"] \n")
        desc = file_describe(filePath,project=project,verbose=False)
        sys.stderr.write(json.dumps(desc,indent=4) + '\n')
    
    return root

def file_find_rep(filePath,project=None,verbose=False,quiet=False):
    '''Returns the replicate tag when in folder {exp_acc}/repN_N.'''

    folder = file_describe(filePath,'folder',project=project,verbose=False)
    rep = ''
    for part in folder.split('/'):
        if part.startswith('rep'):
            rep = part
            break
    if rep == '':
        # Very limited success with file names
        name = file_describe(filePath,'name',verbose=False)
        for part in name.split('_'):
            if part.startswith('rep'):
                rep = part
            elif rep != '':
                if part in ['1','2','3','4','5','6','7','8','9']:
                    rep += '_' + part
                else:
                    break
         
    if verbose:
        sys.stderr.write(rep + '\n')
    if rep == '' and not quiet:  # Note, this is not an error when file is from a different dx context!
        sys.stderr.write("Found nothing for rep: folder["+folder+"] path ["+filePath+"] \n")
        desc = file_describe(filePath,project=project,verbose=False)
        sys.stderr.write(json.dumps(desc,indent=4) + '\n')
    
    return rep

def file_find_exp_id(filePath,project=None,verbose=False,quiet=False):
    '''Returns the experiment id (accession) when in folder {exp_acc}/repN_N.'''

    folder = file_describe(filePath,'folder',project=project,verbose=False)
    exp = ''
    for part in folder.split('/'):
        if part.startswith('ENCSR'):
            exp = part
            break
            
    if exp == '':
        # Very limited success with file names
        name = file_describe(filePath,'name',verbose=False)
        for part in name.split('_'):
            if part.startswith('ENCSR'):
                exp = part
    
    if verbose:
        sys.stderr.write(exp + '\n')
    if exp == '' and not quiet:  # Note, this is not an error when file is from a different dx context!
        sys.stderr.write("Found nothing for exp: folder["+folder+"] path ["+filePath+"] \n")
        desc = file_describe(filePath,project=project,verbose=False)
        sys.stderr.write(json.dumps(desc,indent=4) + '\n')
    
    return exp

def job_describe(job_id,key=None,verbose=False):
    '''Returns dx job's description property matching 'key'.'''

    dxjob = None
    try:
        dxjob = dxpy.get_handler(job_id)
    except:
        sys.stderr.write('ERROR: unable to find job: "' + job_id + '": \n')
        sys.exit(0)  # Do not error on tool run in dx script 
    
    desciption = dxjob.describe()
        
    if not desciption:
        sys.stderr.write('ERROR: unable to find description of job "' + job_id + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script 
    
    if key == None:
        if verbose:
            sys.stderr.write(json.dumps(desciption) + '\n')
        return desciption
    
    if key not in desciption:
        sys.stderr.write('ERROR: unable to find "'+key+'" in description of job "' + job_id + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script
    value = desciption[key]
         
    if verbose:
        sys.stderr.write(value + '\n')
    
    return value

def job_create_root(job_id,verbose=False,quiet=False):
    '''Returns a standard file name root when in folder {exp_acc}/repN_N.'''

    if verbose:
        sys.stderr.write("Trying to create root_name from: "+ job_id + '\n')
    
    folder = job_describe(job_id,'folder',verbose=verbose)
    root = folder_create_root(folder,verbose=verbose)

    if root == '' and not quiet:  # Note, this is not an error when file is from a different dx context!
        sys.stderr.write("Found nothing for root: folder["+folder+"] job ["+job_id+"] \n")
        desc = job_describe(job_id,verbose=False)
        sys.stderr.write(json.dumps(desc,indent=4) + '\n')
    
    return root

def main():
    parser = argparse.ArgumentParser(description =  "Creates a json string of qc_metrics for a given applet. " + \
                                                    "Returns string to stdout and formatted json to stderr.")
    parser.add_argument('-f', '--file',
                        help='DX id, link or path to file.',
                        required=False)
    parser.add_argument('--job',
                        help='DX id of job.',
                        required=False)
    parser.add_argument('-p','--property',
                        help="Property name.",
                        default='QC',
                        required=False)
    parser.add_argument('-s','--subproperty',
                        help="Property name.",
                        default=None,
                        required=False)
    parser.add_argument('-k', '--key',
                        help='Prints just the value for this key.',
                        default=None,
                        required=False)
    parser.add_argument('--keypair',
                        help='Prints the key: value pair for this key.',
                        default=None,
                        required=False)
    parser.add_argument('--project',
                        help="Project (especially helpfule when calling from DX app).",
                        default=None,
                        required=False)
    parser.add_argument('--details',
                        help="Return details json.",
                        default=None, action="store_true", required=False)
    parser.add_argument('--root_name', action="store_true", required=False, default=False, 
                        help="Return a standardized file name root based on file location.")
    parser.add_argument('--rep_tech', action="store_true", required=False, default=False, 
                        help="Return a rep_tech tag based on file location.")
    parser.add_argument('--exp_id', action="store_true", required=False, default=False, 
                        help="Return the exp_id based on file location.")
    parser.add_argument('-d', '--describe', action="store_true", required=False, default=False, 
                        help="Look for key in file description.")
    parser.add_argument('--json', action="store_true", required=False, default=False, 
                        help="Return json.")
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, 
                        help="Suppress non-error stderr messages.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Make some noise.")

    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 2:
        parser.print_usage()
        return
    if args.file == None:
        if args.job == None:
            sys.stderr.write("Requires either '--file' or '--job' argument! \n")
            sys.exit(0)
        elif args.root_name == None:
            sys.stderr.write("Requires '--file' argument! \n")
            sys.exit(0)
    
    
    if args.root_name:
        if args.job:
            root = job_create_root(args.job,verbose=args.verbose,quiet=args.quiet)
        else:
            root = file_create_root(args.file,project=args.project,verbose=args.verbose,quiet=args.quiet)
        print root
        if not args.quiet:
            sys.stderr.write("root_name: '"+root+"'\n")
        sys.exit(0)
        
    elif args.exp_id:
        exp_id = file_find_exp_id(args.file,project=args.project,verbose=args.verbose,quiet=args.quiet)
        print exp_id
        if not args.quiet:
            sys.stderr.write("exp_id: '"+exp_id+"'\n")
        sys.exit(0)
        
    elif args.rep_tech:
        rep = file_find_rep(args.file,project=args.project,verbose=args.verbose,quiet=args.quiet)
        print rep
        if not args.quiet:
            sys.stderr.write("rep: '"+rep+"'\n")
        sys.exit(0)
        
    elif args.details:
        details = file_details(args.file,args.key,project=args.project,verbose=args.verbose)
        print json.dumps(details)
        if args.key != None:
            if not args.quiet:
                sys.stderr.write(args.key + ": '"+json.dumps(details,indent=4,sort_keys=True)+"'\n")
        elif not args.quiet:
            sys.stderr.write(json.dumps(details,indent=4,sort_keys=True)+"\n")
        sys.exit(0)
        
    elif args.describe:
        value = file_describe(args.file,args.key,project=args.project,verbose=args.verbose)
        if args.key != None:
            print value
            if not args.quiet:
                sys.stderr.write(args.key + ": '"+value+"'\n")
            sys.exit(0)
        else:
            properties = value
    else:
        properties = file_get_property(args.file,args.property,args.subproperty,return_json=args.json, \
                                                                        project=args.project,verbose=args.verbose)
        
    # Print out the properties
    if args.key != None:
        if args.key in properties:
            print json.dumps(properties[args.key])
            if not args.quiet:
                sys.stderr.write(json.dumps(properties[args.key],indent=4) + '\n')
        else:
            print ''   
            if not args.quiet:
                sys.stderr.write('(not found)\n')
    elif args.keypair != None:
        if args.keypair in properties:
            print '"' + args.keypair + '": ' + json.dumps(properties[args.keypair])
            if not args.quiet:
                sys.stderr.write('"' + args.keypair + '": ' + json.dumps(properties[args.keypair],indent=4) + '\n')
        else:
            print '"' + args.keypair + '": '
            if not args.quiet:
                sys.stderr.write('"' + args.keypair + '": \n')
    elif isinstance(properties, basestring) or isinstance(properties, int):
        print properties
        if not args.quiet:
            sys.stderr.write(properties + '\n')
    else: 
        print json.dumps(properties)
        if not args.quiet:
            sys.stderr.write(json.dumps(properties,indent=4) + '\n')
    
if __name__ == '__main__':
    main()

