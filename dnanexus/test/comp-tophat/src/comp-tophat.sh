#!/bin/bash
# comp-tophat 0.0.2

#set -x
#set +e

main() {
    #echo "Value of test dataset: '$test_dir'"
    #echo "Value of standard dataset: '$data_dir'"
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/"$test_dir"/*
    #dx download project-BKjyBz00ZZ0PZF5V7Gv00zqG:/test/"$data_dir"/*
    #find

    echo "*****"
    echo "* Running: comp-tophat.sh [aka compare bams]"
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* also using: diff, ls, md5sum, tee, diss and sammy"
    echo "*****"

    set1_bam_fn=`dx describe "$set1_bam" --name`
    set2_bam_fn=`dx describe "$set2_bam" --name`
    bam1=${set1_bam_fn%.bam}
    bam2=${set2_bam_fn%.bam}
    log_diff_fn=compBams_${bam1}_${bam2}.txt
    echo "* Comparing TopHat bam [${bam1}] to [${bam2}]..." | tee ${log_diff_fn}

    echo "* Download files..."
    dx download "$set1_bam"
    dx download "$set2_bam"

    echo " " >>  ${log_diff_fn}
    echo "* Comparing all reads [$set1_bam_fn] to [$set2_bam_fn]..." | tee -a ${log_diff_fn}
    rm -f *.sam
    samtools view -@ 8 ${set1_bam_fn} > ${set1_bam_fn%.bam}.sam 
    samtools view -@ 8 ${set2_bam_fn} > ${set2_bam_fn%.bam}.sam
    #sammy --all ${set1_bam_fn} --fast  # sammy sorts
    #sammy --all ${set2_bam_fn} --fast
    echo "Lines:" | tee -a ${log_diff_fn} 
    wc -l *.sam | tee -a ${log_diff_fn} 
    echo "md5sum:" | tee -a ${log_diff_fn} 
    md5sum *.sam | tee -a ${log_diff_fn}

    echo " " >>  ${log_diff_fn}
    echo "* Comparing uniquely mapped [$set1_bam_fn] to [$set2_bam_fn]..." | tee -a ${log_diff_fn}
    sammy --uniq ${set1_bam_fn} --fast 
    sammy --uniq ${set2_bam_fn} --fast
    echo "Lines:" | tee -a ${log_diff_fn} 
    wc -l *_uniq.sam | tee -a ${log_diff_fn} 
    echo "md5sum:" | tee -a ${log_diff_fn} 
    md5sum *_uniq.sam | tee -a ${log_diff_fn}
    #echo "Split and diff:" | tee -a ${log_diff_fn} 
    #rm -f splitFile?_* 
    #split -l 10000000 ${bam1}_uniq.sam splitFile1_ 
    #split -l 10000000 ${bam2}_uniq.sam splitFile2_ 
    #for f in `ls splitFile1_??`; do
    #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
    #done
    ##diff ${bam1}_uniq.sam ${bam2}_uniq.sam | head | wc -l | tee -a ${log_diff_fn} 


    echo " " >>  ${log_diff_fn}
    echo "* Comparing multi-mapped [$set1_bam_fn] to [$set2_bam_fn]..." | tee -a ${log_diff_fn}
    sammy --multi ${set1_bam_fn} --fast 
    sammy --multi ${set2_bam_fn} --fast 
    echo "Lines:" | tee -a ${log_diff_fn} 
    wc -l *_multi.sam | tee -a ${log_diff_fn} 
    echo "md5sum:" | tee -a ${log_diff_fn} 
    md5sum *_multi.sam | tee -a ${log_diff_fn}
    echo "Split and diff:" | tee -a ${log_diff_fn} 
    rm -f splitFile?_* 
    split -l 10000000 ${bam1}_multi.sam splitFile1_ 
    split -l 10000000 ${bam2}_multi.sam splitFile2_ 
    for f in `ls splitFile1_??`; do
        diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
    done
    ##diff ${bam1}_multi.sam ${bam2}_multi.sam | head | wc -l | tee -a ${log_diff_fn} 


    echo " " >>  ${log_diff_fn}
    echo "* Comparing unmapped [$set1_bam_fn] to [$set2_bam_fn]..." | tee -a ${log_diff_fn}
    sammy --unmapped ${set1_bam_fn} --fast 
    sammy --unmapped ${set2_bam_fn} --fast 
    echo "Lines:" | tee -a ${log_diff_fn} 
    wc -l *_unmapped.sam | tee -a ${log_diff_fn} 
    echo "md5sum:" | tee -a ${log_diff_fn} 
    md5sum *_unmapped.sam | tee -a ${log_diff_fn}
    #echo "Split and diff:" | tee -a ${log_diff_fn} 
    #rm -f splitFile?_* 
    #split -l 10000000 ${bam1}_unmapped.sam splitFile1_ 
    #split -l 10000000 ${bam2}_unmapped.sam splitFile2_ 
    #for f in `ls splitFile1_??`; do
    #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
    #done
    ##diff ${bam1}_unmapped.sam ${bam2}_unmapped.sam | head | wc -l | tee -a ${log_diff_fn} 

    echo " " >>  ${log_diff_fn}
    echo "* Comparing headers [$set1_bam_fn] to [$set2_bam_fn]..." | tee -a ${log_diff_fn}
    sammy --head ${set1_bam_fn} 
    sammy --head ${set2_bam_fn}
    grep -v \@PG ${bam1}_head.sam | grep -v "user command" > ${bam1}_clean_head.sam
    grep -v \@PG ${bam2}_head.sam | grep -v "user command" > ${bam2}_clean_head.sam
    diss ${bam1}_clean_head.sam ${bam2}_clean_head.sam | tee -a ${log_diff_fn}
    #diff ${bam1}_clean_head.sam ${bam2}_clean_head.sam | head | wc -l | tee -a ${log_diff_fn} 
    #diff ${bam1}_clean_head.sam ${bam2}_clean_head.sam | tee -a ${log_diff_fn}


    echo "* Uploading results..."
    log_diff=$(dx upload ${log_diff_fn} --brief)
    dx-jobutil-add-output log_diff "$log_diff" --class=file
    echo "* Finished."
}
