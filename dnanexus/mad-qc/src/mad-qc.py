#!/usr/bin/env python
# Runs "mean absolute deviation" QC metrics on two long-RNA-seq gene quantifications
APP_SCRIPT = "mad-qc.py"
APP_VER = "1.0.1"

import os, subprocess, json
import dxpy

@dxpy.entry_point("main")
def main(quants_a, quants_b):

    # tool_versions.py --applet $script_name --appver $script_ver
    sw_versions = subprocess.check_output(['tool_versions.py', '-a', APP_SCRIPT, '-av', APP_VER])

    dxfile_a = dxpy.DXFile(quants_a)
    dxfile_b = dxpy.DXFile(quants_b)

    print "* Downloading files..."
    dxpy.download_dxfile(dxfile_a.get_id(), "quants_a")
    dxpy.download_dxfile(dxfile_b.get_id(), "quants_b")

    print "* Runnning MAD.R..."
    mad_output = subprocess.check_output(['Rscript', '/usr/bin/MAD.R', 'quants_a', 'quants_b'])
    quants_a_name = dxfile_a.name.split('.')
    quants_b_name = dxfile_b.name.split('.')
    filename = quants_a_name[0] + '_' + quants_b_name[0] + '_' + quants_a_name[1] + '_mad_plot.png'
    subprocess.check_call(['mv', "MAplot.png", filename])
    
    print "* package properties..."
    qc_metrics = {}
    qc_metrics["MAD.R"] = json.loads(mad_output)
    meta_string = json.dumps(qc_metrics)
    print json.dumps(qc_metrics,indent=4)
    props = {}
    props["SW"] = sw_versions

    print "* Upload Plot..."
    plot_dxfile = dxpy.upload_local_file(filename,properties=props,details=qc_metrics)
    
    return { "metadata": meta_string, "mad_plot": plot_dxfile }

dxpy.run()
