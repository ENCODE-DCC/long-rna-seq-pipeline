#!/usr/bin/env python
# Runs "mean absolute deviation" QC metrics on two long-RNA-seq gene quantifications

import os, subprocess, json
import dxpy

def divide_on_common(str_a,str_b):
    '''Divides each string into [common_prefix,variable_middle,common_ending] and returns as set (parts_a,parts_b).'''
    parts_a = ['','','']
    parts_b = ['','','']
    while len(str_a) > 0 and len(str_a) > 0 and str_a[0] == str_b[0]:
        parts_a[0] += str_a[0]
        parts_b[0] += str_b[0]
        str_a = str_a[1:]
        str_b = str_b[1:]
    while len(str_a) > 0 and len(str_a) > 0 and str_a[-1] == str_b[-1]:
        parts_a[2] = str_a[-1] + parts_a[2]
        parts_b[2] = str_b[-1] + parts_b[2]
        str_a = str_a[:-1]
        str_b = str_b[:-1]
    parts_a[1] = str_a
    parts_b[1] = str_b
    return (parts_a,parts_b)

def root_name_from_pair(filename_a,filename_b):
    '''Returns a root name based upon the common and uncommon parts of a pair of file names.'''
    (a_parts, b_parts) = divide_on_common(filename_a,filename_b)
    out_root = ''
    if len(a_parts[0]) > 0:
        out_root = a_parts[0]
    out_root += a_parts[1] + '-' + b_parts[1]
    if len(a_parts[2]) > 0:
        out_root += a_parts[2]
    return out_root # exp1_rep1_1_quants.tsv and exp1_rep2_1_quants.tsv yield exp1_rep1-2_1_quants.tsv
    
@dxpy.entry_point("main")
def main(quants_a, quants_b):

    # tool_versions.py --applet $script_name --appver $script_ver
    sw_versions = subprocess.check_output(['tool_versions.py', '--dxjson', 'dnanexus-executable.json'])

    dxfile_a = dxpy.DXFile(quants_a)
    dxfile_b = dxpy.DXFile(quants_b)

    print "* Downloading files..."
    dxpy.download_dxfile(dxfile_a.get_id(), "quants_a")
    dxpy.download_dxfile(dxfile_b.get_id(), "quants_b")

    # Create and appropriate name for output files
    out_root = root_name_from_pair(dxfile_a.name.split('.')[0],dxfile_b.name.split('.')[0])
    mad_plot_file = out_root + '_mad_plot.png'
        
    # DX/ENCODE independent script is found in resources/usr/bin
    print "* Runnning MAD.R..."
    mad_output = subprocess.check_output(['Rscript', '/usr/bin/MAD.R', 'quants_a', 'quants_b'])
    subprocess.check_call(['mv', "MAplot.png", mad_plot_file ])
    
    print "* package properties..."
    qc_metrics = {}
    qc_metrics["MAD.R"] = json.loads(mad_output)
    meta_string = json.dumps(qc_metrics)
    print json.dumps(qc_metrics,indent=4)
    props = {}
    props["SW"] = sw_versions

    print "* Upload Plot..."
    plot_dxfile = dxpy.upload_local_file(mad_plot_file,properties=props,details=qc_metrics)
    
    return { "metadata": meta_string, "mad_plot": plot_dxfile }

dxpy.run()
