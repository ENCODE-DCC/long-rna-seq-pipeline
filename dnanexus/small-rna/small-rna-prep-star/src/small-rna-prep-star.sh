#!/bin/bash
# small-rna-prep-star.sh

script_name="small-rna-prep-star.sh"
script_ver="2.0.0"

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.2a)
    #(cd STAR; make)

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of ref_genome: '$ref_genome'"
    echo "* Value of annotations: '$annotations'"

    echo "* Download files..."
    ref_root=`dx describe "$ref_genome" --name`
    ref_root=${ref_root%.fasta.gz}
    ref_root=${ref_root%.fa.gz}
    dx download "$ref_genome" -o "$ref_root".fa.gz
    gunzip "$ref_root".fa.gz

    annotation_root=`dx describe "$annotations" --name`
    annotation_root=${annotation_root%.gtf.gz}
    dx download "$annotations" -o "$annotation_root".gtf.gz
    gunzip ${annotation_root}.gtf.gz

    echo "* Reference file(s): '${ref_root}.fa'"
    
    # hg19/mm10 and male/female?
    if [ -f /usr/bin/parse_property.py ]; then
        genome=`parse_property.py -f "$ref_genome" -p genome`
        gender=`parse_property.py -f "$ref_genome" -p gender`
        anno=`parse_property.py -f "$annotations" -p annotation`
    fi
    if [ "$genome" == "" ]; then
        if [[ $genome_root == *"hg19"* ]]; then
            genome="hg19"
        elif [[ $genome_root == *"hg38"* ]]; then
            genome="hg38"
        elif [[ $genome_root == *"mm10"* ]]; then
            genome="mm10"
        else
            genome="unknown"
        fi
    fi
    if [ "$gender" == "" ]; then
        if [[ $genome_root == *"female"* ]]; then
            gender="female"
        elif [[ $genome_root == *"male"* ]]; then
            gender="male"
        else
            gender="unknown"
        fi
    fi
    if [ "$anno" == "" ]; then
        if [[ $annotation_root == *"v19"* ]] || [[ $annotation_root == *"V19"* ]]; then
            anno="v19"
        elif [[ $annotation_root == *"M4"* ]]; then
            anno="M4"
        elif [[ $annotation_root == *"M3"* ]]; then
            anno="M3"
        elif [[ $annotation_root == *"M2"* ]]; then
            anno="M2"
        else
            anno="unknown"
        fi
    fi
    archive_root="${genome}_${gender}_${anno}_sRNA_starIndex"
    
    echo "* Build index for '${genome}-${gender} and annotation ${anno}'..."
    set -x
    mkdir out
    STAR --runMode genomeGenerate --genomeFastaFiles ${ref_root}.fa --sjdbGTFfile ${annotation_root}.gtf \
            --sjdbOverhang 1 --runThreadN 8 --genomeDir out/ --outFileNamePrefix out
    set +x
    
    echo "* Build gene ID file from the annotation for later use..."
    # Create a geneID file and put it in the index archive for later use:
    set -x
    awk '$3=="gene" || substr($14,2,length($14)-3)=="tRNAscan" {g=substr($14,2,length($14)-3); if (g=="miRNA" || g=="snoRNA" || g=="snRNA" || g=="tRNAscan") {print substr($10,2,length($10)-3)} }' \
        ${annotation_root}.gtf | sort > out/smallRNA.geneID
    set +x

    # Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    set -x
    refComment="@CO\tREFID:$(basename ${ref_root})"
    echo -e ${refComment} > out/star_bamCommentLines.txt
    echo `cat "out/star_bamCommentLines.txt"`
    set +x

    echo "* Tar and upload results..."
    echo "out"
    set -x
    echo `ls out/`
    tar -czf ${archive_root}.tgz out/
    set +x

    star_index=$(dx upload ${archive_root}.tgz --property genome="$genome" --property gender="$gender" \
                                               --property annotation="$anno" --property SW="$versions" --brief)
    dx-jobutil-add-output star_index $star_index --class=file

    echo "* Finished."
}
