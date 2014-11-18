#!/bin/bash
# rampage-peaks 0.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0d)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: rampage-peaks.sh [v0.0.1]"
    echo "* STAR version:     ["`STAR --version | awk '{print $1}' | cut -d _ -f 2-`"]"
    #echo "* bedGraphToBigWig version: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "* wigToBigWig version: "`wigToBigWig 2>&1 | grep "wigToBigWig v" | awk '{print $2$3}'`
    echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{printf "v%s", $3}'`
    # TODO: rampagePeakCaller.py version?
    echo "*****"

    echo "Value of marked_bam: '$marked_bam'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_fn=`dx describe "$marked_bam" --name`
    bam_fn=${bam_fn%.bam}
    echo "* Bam file: '"$bam_fn".bam'"
    dx download "$marked_bam" -o "$bam_fn".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    echo "* Making 5' signal..."
    mkdir -p Signal
    STAR --inputBAMfile ${bam_fn}.bam --runMode inputAlignmentsFromBAM \
         --outWigReferencesPrefix chr --outWigType wiggle read1_5p --outWigNorm None \
         --outFileNamePrefix Signal/read1_5p.
    #mv Signal/Signal*bg .

    echo "* Concatenating 5p strands, adding track lines..."
    echo track type=wiggle_0 name=read1_5p_+ > read1_5p.wig
    cat Signal/read1_5p.Signal.Unique.str1.out.wig >> read1_5p.wig
    echo track type=wiggle_0 name=read1_5p_- >> read1_5p.wig
    cat Signal/read1_5p.Signal.Unique.str2.out.wig >> read1_5p.wig

    echo "* Making downstream read signal..."
    STAR --inputBAMfile ${bam_fn}.bam --runMode inputAlignmentsFromBAM \
         --outWigReferencesPrefix chr --outWigType wiggle read2 --outWigNorm None \
         --outFileNamePrefix Signal/read2.

    echo "* Concatenating downstream reads, adding track lines..."
    echo track type=wiggle_0 name=read2_+ > read2.wig
    cat Signal/read2.Signal.Unique.str1.out.wig >> read2.wig
    echo track type=wiggle_0 name=read2_- >> read2.wig
    cat Signal/read2.Signal.Unique.str2.out.wig >> read2.wig

    echo "* Converting wigs to bigWigs..."
    # TODO: What signal files are really being kept?
    #bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_fn}_minusAll.bw
    #bedGraphToBigWig Signal.Unique.str1.out.bg         chromSizes.txt ${bam_fn}_minusUniq.bw
    #bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg chromSizes.txt ${bam_fn}_plusAll.bw
    #bedGraphToBigWig Signal.Unique.str2.out.bg         chromSizes.txt ${bam_fn}_plusUniq.bw
    wigToBigWig read1_5p.wig chromSizes.txt ${bam_fn}_read1_5p.bw
    wigToBigWig read2.wig chromSizes.txt ${bam_fn}_read2.bw
    echo `ls`

    echo "* Calling peaks..."
    python2.7 /usr/bin/rampagePeakCaller.py --threads 8 --exclude chrM \
        --read2-background read2.wig \
        --bg-weight 0.1  --NBinom 0.6 --fdr 0.00000001 chrNL.txt read1_5p.wig
 
    echo "* Canverting bed 6 to bigWig..."
    # TODO what is the name of the rampagePeakCaller.py created bed6 file?
    mv rampagePeaks.bed ${bam_fn}_rampage_peaks.bed
    bedToBigBed ${bam_fn}_rampage_peaks.bed chromSizes.txt ${bam_fn}_rampage_peaks.bb

    echo "* Upload results..."
    reads1_5p_bw=$(dx upload ${bam_fn}_read1_5p.bw --brief)
    reads2_bw=$(dx upload ${bam_fn}_read2.bw --brief)
    rampage_peaks=$(dx upload ${bam_fn}_rampage_peaks.bed --brief)
    rampage_peaks_bb=$(dx upload ${bam_fn}_rampage_peaks.bb --brief)

    dx-jobutil-add-output reads1_5p_bw "$reads1_5p_bw" --class=file
    dx-jobutil-add-output reads2_bw "$reads2_bw" --class=file
    dx-jobutil-add-output rampage_peaks "$rampage_peaks" --class=file
    dx-jobutil-add-output rampage_peaks_bb "$rampage_peaks_bb" --class=file

    #echo "* Temporary uploads..."
    # temprary for comparison only!
    #mv Signal.UniqueMultiple.str1.out.bg ${bam_fn}_minusAll.bg
    #mv Signal.Unique.str1.out.bg         ${bam_fn}_minusUniq.bg
    #mv Signal.UniqueMultiple.str2.out.bg ${bam_fn}_plusAll.bg
    #mv Signal.Unique.str2.out.bg         ${bam_fn}_plusUniq.bg
    #minusAll_bg=$(dx upload ${bam_fn}_minusAll.bg --brief)
    #plusAll_bg=$(dx upload ${bam_fn}_plusAll.bg --brief)
    #minusUniq_bg=$(dx upload ${bam_fn}_minusUniq.bg --brief)
    #plusUniq_bg=$(dx upload ${bam_fn}_plusUniq.bg --brief)
    #dx-jobutil-add-output minusAll_bg "$minusAll_bg" --class=file
    #dx-jobutil-add-output plusAll_bg "$plusAll_bg" --class=file
    #dx-jobutil-add-output minusUniq_bg "$minusUniq_bg" --class=file
    #dx-jobutil-add-output plusUniq_bg "$plusUniq_bg" --class=file
    echo "* Finished."
}
