#!/bin/bash
# comp-signals-se 0.0.2

#set -x
#set +e

main() {
    echo "*****"
    echo "* Running: comp-signals-se.sh"
    echo "* Using: awk, tee, alexAwkBg.sh"
    echo "*****"

    # TODO: at some point the bam-to-bw scripts will stop uploading the bedGraphs.
    #       When that happens, this compare routine can be switched to take bigWigs and
    #       use bigWigToBedGraph before the comparisons.

    set1_all_fn=`dx describe "$set1_all" --name`
    set2_all_fn=`dx describe "$set2_all" --name`
    set1=${set1_all_fn%_all.bg}
    set2=${set1_all_fn%_all.bg}
    log_diff_fn=compBg_${set1}_${set2}.txt
    echo "* Comparing bedGraphs [${set1}] to [${set2}]..." | tee ${log_diff_fn}
    echo " " >>  ${log_diff_fn}

    echo "* Comparing BG All [$set1_all_fn] to [$set2_all_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_all"
    dx download "$set2_all"
    alexAwkBg.sh ${set1_all_fn} ${set2_all_fn} | grep Maximum | tee -a ${log_diff_fn} 


    echo " " >>  ${log_diff_fn}
    set1_uniq_fn=`dx describe "$set1_uniq" --name`
    set2_uniq_fn=`dx describe "$set2_uniq" --name`
    dx download "$set1_uniq"
    dx download "$set2_uniq"
    echo "* Comparing BG Unique [$set1_uniq_fn] to [$set2_uniq_fn]..." | tee -a ${log_diff_fn}
    alexAwkBg.sh ${set1_uniq_fn} ${set2_uniq_fn} | grep Maximum | tee -a ${log_diff_fn} 


    #echo "Value of test dataset: '$test_dir'"
    #echo "Value of standard dataset: '$data_dir'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe"
    # "$variable" --name".

    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/"$test_dir"/*
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/test/"$data_dir"/*
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:test/"$data_dir"/Log.final.out
    #find

    echo "* Uploading results..."
    log_diff=$(dx upload ${log_diff_fn} --brief)
    dx-jobutil-add-output log_diff "$log_diff" --class=file
    echo "* Finished."
}
