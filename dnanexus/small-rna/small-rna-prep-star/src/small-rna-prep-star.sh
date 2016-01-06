#!/bin/bash
# small-rna-prep-star.sh

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.2a)
    #(cd STAR; make)

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of ref_genome:  '$ref_genome'"
    echo "* Value of annotations: '$annotations'"

    echo "* Download files..."
    ref_root=`dx describe "$ref_genome" --name`
    ref_root=${ref_root%.fasta.gz}
    ref_root=${ref_root%.fa.gz}
    dx download "$ref_genome" -o ${ref_root}.fa.gz

    anno_root=`dx describe "$annotations" --name`
    anno_root=${anno_root%.gtf.gz}
    dx download "$annotations" -o ${anno_root}.gtf.gz

    echo "* Reference file: '${ref_root}.fa.gz'"
    echo "* Annotation file: '${anno_root}.gtf.gz'"
    
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
    # Make certain the extraction script is in place
    cp /usr/bin/extract_gene_ids.awk .
    
    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    srna_index.sh ${ref_root}.fa.gz ${anno_root}.gtf.gz $anno $genome $gender
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    archive_root="${genome}_${anno}_sRNA_starIndex"
    if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
        archive_root="${genome}_${gender}_${anno}_sRNA_starIndex"
    fi
    gene_id_root="${genome}_${anno}_gene_ids"

    echo "Upload results..."
    star_index=$(dx upload ${archive_root}.tgz --property genome="$genome"   --property gender="$gender" \
                                               --property annotation="$anno" --property SW="$versions" --brief)
    star_genes=$(dx upload ${gene_id_root}.txt --property genome="$genome"  --property annotation="$anno" \
                                                                            --property SW="$versions" --brief)
    dx-jobutil-add-output star_index $star_index --class=file
    dx-jobutil-add-output star_genes $star_genes --class=file

    echo "* Finished."
}
