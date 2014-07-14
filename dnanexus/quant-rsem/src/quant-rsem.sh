#!/bin/bash
# quant-rsem 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {

    echo "Value of annotation_bam: '$annotation_bam'"
    echo "Value of read_prefix: '$read_prefix'"
    echo "Value of paired_end: '$paired_end'"
    echo "Value of stranded: '$stranded'"
    echo "Value of rsem_index: '$rsem_index'"
    echo "Random number seed": '$rnd_seed'""

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    dx download "$annotation_bam" -o annotation.bam

    dx download "$rsem_index" -o rsem_index.tgz
    tar zxvf rsem_index.tgz

    index_prefix=`ls out/*.grp | cut -d'.' -f1`
    echo "found $index_prefix grp file"

    echo "dowload and install RSEM"
    git clone https://github.com/bli25wisc/RSEM.git
    git checkout 92b24279a3ecc72946e7e7c23149ad0d181f373a
    ## get correct commit bit from submodule
    (cd RSEM; make)

    extraFlags=""
    if [ -n "$paired" ]
    then
        extraFlags="--paired-end "
    fi

    if [ -n "stranded" ]
    then
        extraFlags=${extraFlags}"--forward-prob 0"
    fi

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.
    RSEM/rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} \
                                 -p 8 --ci-memory 60000 ${extraFlags} \
                                 annotation.bam ${index_prefix} ${read_prefix}_rsem_quant

    echo `ls ${read_prefix}*`
    # deliver results:
    #mv rsemOut.genes.results ${read_prefix}.genes.rsem.results
    #mv rsemOut.isoforms.results ${read_prefix}.isoforms.rsem.results
    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    genomic_quant=$(dx upload ${read_prefix}_rsem_quant.genes.results --brief)
    transcript_quant=$(dx upload ${read_prefix}_rsem_quant.isoforms.results --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output genomic_quant "$genomic_quant" --class=file
    dx-jobutil-add-output transcript_quant "$transcript_quant" --class=file

}
