#!/bin/bash
# prep-rsem.sh

main() {
    # Now in resources/usr/bin
    #echo "* Dowload and install RSEM..."
    #git clone https://github.com/bli25wisc/RSEM.git
    # (cd RSEM; git checkout tags/v1.2.19)
    ### get correct commit bit from submodule
    #(cd RSEM; make)

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of ref_genome:  '$ref_genome'"
    echo "* Value of spike_in:    '$spike_in'"
    echo "* Value of annotations: '$annotations'"

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
    lrna_index_rsem.sh ${ref_root}.fa.gz ${spike_root}.fa.gz ${anno_root}.gtf.gz $anno $genome $gender
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    archive_file="${genome}_${anno}_${spike_root}_rsemIndex.tgz"
    if [ "$gender" == "female" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
        archive_file="${genome}_${gender}_${anno}_${spike_root}_rsemIndex.tgz"
    fi

    echo "* Upload results..."
    rsem_index=$(dx upload $archive_file --property genome="$genome"   --property gender="$gender" \
                                         --property annotation="$anno" --property spike_in="$spike_root" \
                                         --property SW="$versions" --brief)

    dx-jobutil-add-output rsem_index $rsem_index --class=file
    echo "* Finished."
}
