#!/bin/bash
# prep-tophat.sh

main() {
    # Now in resources/usr/bin
    ## install tophat 2.0.8
    #wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.8.Linux_x86_64.tar.gz
    ## install bowtie2_2.1.0
    #wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of ref_genome:  '$ref_genome'"
    echo "* Value of spike_in:    '$spike_in'"
    echo "* Value of annotations: '$annotations'"
    echo "* Value of tiny_fq:     '$tiny_fq'"

    echo "* Value of annotations: '$annotations'"
    echo "* Value of genome:   '$genome'"

    echo "* Download files..."
    echo "* Download files..."
    ref_root=`dx describe "$ref_genome" --name`
    ref_root=${ref_root%.fasta.gz}
    ref_root=${ref_root%.fa.gz}
    dx download "$ref_genome" -o ${ref_root}.fa.gz

    spike_root=`dx describe "$spike_in" --name`
    spike_root=${spike_root%.fasta.gz}
    spike_root=${spike_root%.fa.gz}
    dx download "$spike_in" -o ${spike_root}.fa.gz

    anno_root=`dx describe "$annotations" --name`
    anno_root=${anno_root%.gtf.gz}
    dx download "$annotations" -o ${anno_root}.gtf.gz


    dx download "$tiny_fq" -o tiny.fq.gz

    echo "* Reference file(s): '$ref'"
    echo "* Value of geno_prefix: '$geno_prefix'"
    echo "* Value of anno_prefix: '$anno_prefix'"

    # hg19/mm10 and male/female?
    if [ -f /usr/bin/parse_property.py ]; then
        genome=`parse_property.py -f "$ref_genome" -p genome`
        gender=`parse_property.py -f "$ref_genome" -p gender`
        anno=`parse_property.py -f "$annotations" -p annotation`
    fi
    if [ "$genome" == "" ]; then
        if [[ $ref_root == *"hg19"* ]]; then
            genome="hg19"
        elif [[ $ref_root == *"GRCh38"* ]]; then
            genome="GRCh38"
        elif [[ $ref_root == *"mm10"* ]]; then
            genome="mm10"
        fi
    fi
    if [ "$genome" == "" ]; then
        genome="unknown"
        echo "* WARNING genome: '$genome'" 
    else
        echo "* genome: '$genome'" 
    fi
    if [ "$gender" == "" ]; then
        if [[ $ref_root == *"female"* ]] || [[ $ref_root == *"XX"* ]]; then
            gender="XX"
        elif [[ $ref_root == *"male"* ]] || [[ $ref_root == *"XY"* ]]; then
            gender="XY"
        fi
    fi
    if [ "$gender" != "" ]; then
        echo "* gender: '$gender'" 
    fi
    if [ "$anno" == "" ]; then
        if [[ $anno_root == *"v19"* ]] || [[ $anno_root == *"V19"* ]]; then
            anno="v19"
            echo "* WARNING annotation version: '$anno'" 
        elif [[ $anno_root == *"v24"* ]] || [[ $anno_root == *"V24"* ]]; then
            anno="v24"
        elif [[ $anno_root == *"M4"* ]]; then
            anno="M4"
        elif [[ $anno_root == *"M3"* ]]; then
            anno="M3"
        elif [[ $anno_root == *"M2"* ]]; then
            anno="M2"
        fi
    fi
    if [ "$anno" == "" ]; then
        anno="unknown"
        echo "* WARNING annotation version: '$anno'"
    else 
        echo "* annotation version: '$anno'" 
    fi

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    lrna_index_tophat.sh ${ref_root}.fa.gz ${spike_root}.fa.gz ${anno_root}.gtf.gz tiny.fq.gz $anno $genome $gender
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    archive_file="${genome}_${anno}_${spike_root}_tophatIndex.tgz"
    if [ "$gender" == "female" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
        archive_file="${genome}_${gender}_${anno}_${spike_root}_tophatIndex.tgz"
    fi

    echo "* Upload results..."
    tophat_index=$(dx upload $archive_file --property genome="$genome"   --property gender="$gender" \
                                           --property annotation="$anno"  --property spike_in="$spike_root" \
                                           --property SW="$versions" --brief)

    dx-jobutil-add-output tophat_index $tophat_index --class=file
    echo "* Finished."
}
