#!/bin/bash
# align-tophat-pe.sh

main() {
    # Now in resources/usr/bin
    ## install tophat 2.0.8
    #wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.8.Linux_x86_64.tar.gz
    ## install bowtie2_2.1.0
    #wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz
    #wget wget https://github.com/xweigit/xweiEncodeScripts/archive/v1.0.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
    #echo "*****"
    #echo "* Running: align-tophat-pe.sh v[1.0.1]"
    #echo "* TopHat version: "`tophat -v | awk '{print $2}'`
    #echo "* bowtie2 version: "`bowtie2 --version 2>&1 | grep bowtie | awk '{print $3}'`
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    #echo "* tophat_bam_xsA_tag_fix.pl version: "`perl /usr/bin/tophat_bam_xsA_tag_fix.pl --version 2>&1`
    #echo "*****"

    echo "* Value of reads: '$reads_1'"
    echo "* Value of reads: '$reads_2'"
    echo "* Value of tophat_index: '$tophat_index'"
    echo "* Value of library_id: '$library_id'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads1[@]}
    do
        file_root=`dx describe "${reads1[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${concat}" == "" ]; then
                outfile_name="${outfile_name}_concat" 
                concat="s concatenated as"
            fi
        fi
        echo "* Downloading concatenating ${file_root}.fq.gz file..."
        dx download "${reads1[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Reads1 fastq${concat} file: '${outfile_name}.fq.gz'"
    reads1_root=${outfile_name}
    ls -l ${reads1_root}.fq.gz

    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads2[@]}
    do
        file_root=`dx describe "${reads2[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${concat}" == "" ]; then
                outfile_name="${outfile_name}_concat" 
                concat="s concatenated as"
            fi
        fi
        echo "* Downloading and concatenating ${file_root}.fq.gz file..."
        dx download "${reads2[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Reads2 fastq${concat} file: '${outfile_name}.fq.gz'"
    ls -l ${outfile_name}.fq.gz
    reads2_root=${outfile_name}
    ls -l ${reads2_root}.fq.gz
    bam_root="${reads1_root}_${reads2_root}_tophat"
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py -f "'${reads1[0]}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            bam_root="${new_root}_tophat"
        fi
    fi
    echo "* Alignments file will be: '${bam_root}.bam'"

    echo "* Downloading and extracting TopHat index archive..."
    dx download "$tophat_index" -o tophat_index.tgz
    tar zxvf tophat_index.tgz

    # unzips into "out/"
    gff=`ls out/*.gff`
    anno_prefix=${gff%.gff}
    bamComments=`ls out/*_bamCommentLines.txt`
    geno_prefix=${bamComments%_bamCommentLines.txt}
    echo "* Value of geno_prefix: '$geno_prefix'"
    echo "* Value of anno_prefix: '$anno_prefix'"

    # Fill in your application code here.

    echo "* Map reads..."
    set -x
    tophat -p ${nthreads} -z0 -a 8 -m 0 --min-intron-length 20 --max-intron-length 1000000 \
        --read-edit-dist 4 --read-mismatches 4 -g 20  --no-discordant --no-mixed \
        --library-type fr-firststrand --transcriptome-index ${anno_prefix} \
        ${geno_prefix} ${reads1_root}.fq.gz ${reads2_root}.fq.gz
    set +x
    ls -l tophat_out/accepted_hits.bam
    ls -l tophat_out/unmapped.bam

    echo "* Set up headers..."
    set -x
    HD="@HD\tVN:1.4\tSO:coordinate" 
    stCommand="perl tophat_bam_xsA_tag_fix.pl tophat_out/accepted_hits.bam | samtools view -bS - | samtools sort - mapped_fixed; samtools merge -h newHeader.sam merged.bam mapped_fixed.bam out/unmapped.bam"
    newPG="@PG\tID:Samtools\tPN:Samtools\tCL:"$stCommand"\tPP:Tophat\tVN:VN:0.1.17 (r973:277)"
    libraryComment="@CO\tLIBID:${library_id}"

    samtools view -H tophat_out/accepted_hits.bam | \
    gawk -v HD="$HD" -v newPG="$newPG" -v library="$libraryComment" \
        '{     if ($0 ~ /^@PG/) {PG=$0} 
          else if ($0 ~ /^@HD/) {print HD; }
          else if($0 ~ /^@SQ/) {print $0};
         }; END{print newPG"\n"PG"\n"library;}' > newHeader.sam

    # Add reference genome and transcriptome used
    cat ${geno_prefix}_bamCommentLines.txt >> newHeader.sam
    set +x

    echo "* Fix unmapped bam and sort before merge..."
    set -x
    perl /usr/bin/tophat_bam_xsA_tag_fix.pl tophat_out/accepted_hits.bam | \
                samtools view -bS - | samtools sort - mapped_fixed
    set +x
    ls -l mapped_fixed.bam

    echo "* Merge aligned and unaligned into single bam, using the patched up header..."
    set -x
    samtools merge -h newHeader.sam merged.bam mapped_fixed.bam tophat_out/unmapped.bam
    mv merged.bam ${bam_root}.bam
    set +x
    ls -l ${bam_root}.bam

    echo "* Upload results..."
    tophat_bam=$(dx upload ${bam_root}.bam --property SW="$versions" --brief)
    dx-jobutil-add-output tophat_bam "$tophat_bam" --class=file
    echo "* Finished."
}

