#!/bin/bash
# comp-star-rsem 0.0.2

#set -x
#set +e

main() {
    echo "*****"
    echo "* Running: comp-star-rsem.sh"
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* also using: cut, diff, ls, md5sum, tee, diss"
    echo "*****"

    set1_log_final_fn=`dx describe "$set1_log_final" --name`
    set2_log_final_fn=`dx describe "$set2_log_final" --name`
    set1=${set1_log_final_fn%_Log.final.out}
    set2=${set1_log_final_fn%_Log.final.out}
    log_diff_fn=compStarRsem_${set1}_${set2}.txt
    echo "* Comparing STAR/RSEM [${set1}] to [${set2}]..." | tee ${log_diff_fn}

    echo " " >>  ${log_diff_fn}
    echo "* Comparing STAR [$set1_log_final_fn] to [$set2_log_final_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_log_final"
    dx download "$set2_log_final"
    diss <(awk 'NR>4{print}' ${set1_log_final_fn}) <(awk 'NR>4{print}' ${set2_log_final_fn}) | tee -a ${log_diff_fn}


    echo " " >>  ${log_diff_fn}
    set1_genome_bam_fn=`dx describe "$set1_genome_bam" --name`
    set2_genome_bam_fn=`dx describe "$set2_genome_bam" --name`
    echo "* Comparing STAR Genome aligned bam [$set1_genome_bam_fn] to [$set2_genome_bam_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_genome_bam"
    dx download "$set2_genome_bam"
    samtools view -@ 8 ${set1_genome_bam_fn} > ${set1_genome_bam_fn%.bam}_geno.sam 
    samtools view -@ 8 ${set2_genome_bam_fn} > ${set2_genome_bam_fn%.bam}_geno.sam
    echo "- Lines:" | tee -a ${log_diff_fn} 
    wc -l *_geno.sam | tee -a ${log_diff_fn} 
    echo "- md5sum:" | tee -a ${log_diff_fn} 
    md5sum *_geno.sam | tee -a ${log_diff_fn}
    #echo "Split and diff:" | tee -a ${log_diff_fn}
    #rm -f splitFile?_* 
    #split -l 10000000 ${set1_genome_bam_fn%.bam}_geno.sam splitFile1_ 
    #split -l 10000000 ${set2_genome_bam_fn%.bam}_geno.sam splitFile2_ 
    #for f in `ls splitFile1_??`; do
    #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
    #done
    ##diss <(samtools view ${set1_genome_bam_fn}) <(samtools view ${set2_genome_bam_fn}) | tee -a ${log_diff_fn}


    echo " " >>  ${log_diff_fn}
    set1_anno_bam_fn=`dx describe "$set1_anno_bam" --name`
    set2_anno_bam_fn=`dx describe "$set2_anno_bam" --name`
    echo "* Comparing STAR Annotation aligned bam [$set1_anno_bam_fn] to [$set2_anno_bam_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_anno_bam"
    dx download "$set2_anno_bam"
    samtools view -@ 8 ${set1_anno_bam_fn} > ${set1_anno_bam_fn%.bam}_anno.sam 
    samtools view -@ 8 ${set2_anno_bam_fn} > ${set2_anno_bam_fn%.bam}_anno.sam
    echo "- Lines:" | tee -a ${log_diff_fn} 
    wc -l *_anno.sam | tee -a ${log_diff_fn} 
    echo "- md5sum:" | tee -a ${log_diff_fn} 
    md5sum *_anno.sam | tee -a ${log_diff_fn}
    #echo "Split and diff:" | tee -a ${log_diff_fn} 
    #rm -f splitFile?_* 
    #split -l 10000000 ${set1_anno_bam_fn%.bam}_anno.sam splitFile1_ 
    #split -l 10000000 ${set2_anno_bam_fn%.bam}_anno.sam splitFile2_ 
    #for f in `ls splitFile1_??`; do
    #    diss $f splitFile2_${f#splitFile1_} | tee -a ${log_diff_fn}
    #done
    ##diss <(samtools view ${set1_anno_bam_fn}) <(samtools view ${set2_anno_bam_fn}) | tee -a ${log_diff_fn}


    echo " " >>  ${log_diff_fn}
    set1_isoform_results_fn=`dx describe "$set1_isoform_results" --name`
    set2_isoform_results_fn=`dx describe "$set2_isoform_results" --name`
    echo "* Comparing RSEM isoforms results [$set1_isoform_results_fn] to [$set2_isoform_results_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_isoform_results"
    dx download "$set2_isoform_results"
    #cut -f1-8 Quant.isoforms.results > iso.a.diff
    #cut -f1-8 rsem_isoform_quant > iso.b.diff
    #echo `ls *diff`
    #diff iso.a.diff iso.b.diff > isoform_quant_diff
    diss <(cut -f1-8 ${set1_isoform_results_fn}) <(cut -f1-8 ${set2_isoform_results_fn}) | tee -a ${log_diff_fn}


    echo " " >>  ${log_diff_fn}
    set1_gene_results_fn=`dx describe "$set1_gene_results" --name`
    set2_gene_results_fn=`dx describe "$set2_gene_results" --name`
    echo "* Comparing RSEM genes results [$set1_gene_results_fn] to [$set2_gene_results_fn]..." | tee -a ${log_diff_fn}
    dx download "$set1_gene_results"
    dx download "$set2_gene_results"
    #for ii in `cd $data_dir; ls *bw`
    #do
    #    echo $ii
    #    diff $data_dir/$ii $2/$ii | head
    #done
    diss <(cut -f1-7 ${set1_gene_results_fn}) <(cut -f1-7 ${set2_gene_results_fn}) | tee -a ${log_diff_fn}

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
