#!/bin/bash
# rampage-idr 0.1.0

main() {
    echo "* Installing Anacaonda3 (python3.4.3, numpy-1.9.2 matplotlib-1.4.3..."
    set -x
    wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda3-2.2.0-Linux-x86_64.sh >> ../install.log 2>&1
    bash Anaconda3-2.2.0-Linux-x86_64.sh -b
    ana_bin=`echo ~/anaconda3/bin`
    # python symlink will interfere with python2.7
    rm ${ana_bin}/python
    ls ${ana_bin}
    export PATH=${ana_bin}:$PATH
    python3 -V 2>&1 | tee -a install.log
    #sudo python3 -V 2>&1 | tee -a install.log
    #sudo ${ana_bin}/python3 -V 2>&1 | tee -a install.log
    set +x

    echo "* Installing idr..."
    set -x
    wget https://github.com/nboley/idr/archive/2.0.0beta5.tar.gz -O idr.tgz >> ../install.log 2>&1
    mkdir idr
    tar -xzf idr.tgz -C idr --strip-components=1
    cd idr
    #sudo python3 setup.py install >> ../install.log 2>&1
    sudo ${ana_bin}/python3 setup.py install
    cd ..
    set +x
    
    echo "*****"
    echo "* Running: rampage-idr.sh [v0.1.0]"
    echo "* idr version: "`idr/bin/idr --version 2>&1 | grep IDR | awk '{print $2}'`
    echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{printf "v%s", $3}'`
    echo "*****"

    echo "Value of peaks_a:     '$peaks_a'"
    echo "Value of peaks_b:     '$peaks_b'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    peaks_a_fn=`dx describe "$peaks_a" --name`
    peaks_a_fn=${peaks_a_fn%_rampage_peaks.bed}
    peaks_a_fn=${peaks_a_fn%.bed}
    dx download "$peaks_a" -o peaks_a.bed
    echo "* First bed: '"$peaks_a_fn".bed'"

    peaks_b_fn=`dx describe "$peaks_b" --name`
    peaks_b_fn=${peaks_b_fn%_rampage_peaks.bed}
    peaks_b_fn=${peaks_b_fn%.bed}
    dx download "$peaks_b" -o peaks_b.bed
    echo "* Second bed: '"$peaks_b_fn".bed'"

    dx download "$chrom_sizes" -o chromSizes.txt

    idr_root=${peaks_a_fn}_${peaks_b_fn}_idr
    echo "* Rampage IDR root: '"$idr_root"'"

    # TODO:
    # Proper paramters for idr
    # Should we remove the ERCC spike-in peaks first?
    # Do we want bed?  bigBed?  Limit to passing peaks only? (--idr-threshold 0.05)
    # Do we want plot?
    
    echo "* Removing any spike-ins from bed files..."
    set -x
    grep "^chr" peaks_a.bed > peaks_a_clean.bed
    grep "^chr" peaks_b.bed > peaks_b_clean.bed
    set -x

    echo "* Running IDR..."
    # Note: --plot: "WARNING: matplotlib does not appear to be installed and is required for plotting"
    #               even though "sudo apt-get install python-matplotlib" says it is installed 
    set -x
    idr/bin/idr --input-file-type bed --rank 7 --plot --samples peaks_a_clean.bed peaks_b_clean.bed 2>&1 | tee idr_summary.txt
    sort -k1,1 -k2,2n < idrValues.txt > ${idr_root}.bed
    touch idrValues.txt.png
    mv idrValues.txt.png ${idr_root}.png
    set +x

    echo "* Converting bed to bigBed..."
    set -x
    bedToBigBed ${idr_root}.bed -type=bed6+ -as=/usr/bin/rampage_idr_peaks.as chromSizes.txt ${idr_root}.bb
    set +x

    # Gather metrics
    meta=`echo \"IDR summary\": { `
    #          Initial parameter values: [0.10 1.00 0.20 0.50]
    var=`grep "Initial parameter values" idr_summary.txt | awk '{printf "%s, %s, %s, %s",$4,$5,$6,$7}'`
    var=`echo \"Initial mu sigma rho and mix values\": $var`
    meta=`echo $meta $var`
    #          Final parameter values: [0.09 0.20 0.10 0.99]
    var=`grep "Final parameter values" idr_summary.txt | awk '{printf "%s, %s, %s, %s",$4,$5,$6,$7}'`
    var=`echo \"Final mu sigma rho and mix values\": $var`
    meta=`echo $meta, $var`
    #          Number of reported peaks - 53/53 (100.0%)
    var=`grep "Number of reported peaks" idr_summary.txt | awk '{print $6}' | cut -d / -f 2`
    var=`echo \"Total number of peaks\": $var`
    meta=`echo $meta, $var`
    #          Number of peaks passing IDR cutoff of 0.05 - 41/53 (77.4%)
    var=`grep "Number of peaks passing IDR cutoff" idr_summary.txt | awk '{print $8}'`
    var=`echo \"IDR cutoff\": $var`
    meta=`echo $meta, $var`
    var=`grep "Number of peaks passing IDR cutoff" idr_summary.txt | awk '{print $10}' | cut -d / -f 1`
    var=`echo \"Number of peaks passing IDR cutoff\": $var`
    meta=`echo $meta, $var`
    var=`grep "Number of peaks passing IDR cutoff" idr_summary.txt | awk '{print $11}' | tr -d \(\)\% | awk '{print $1}'`
    var=`echo \"Percent of peaks passing IDR cutoff\": $var`
    meta=`echo $meta, $var }`
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    details=`echo { $meta }`
    rampage_idr_bed=$(dx upload ${idr_root}.bed --details="$details" --property QC="$meta" --brief)
    rampage_idr_bb=$(dx upload ${idr_root}.bb --details="$details" --property QC="$meta" --brief)
    rampage_idr_png=$(dx upload ${idr_root}.png --details="$details" --property QC="$meta" --brief)

    dx-jobutil-add-output rampage_idr_bed "$rampage_idr_bed" --class=file
    dx-jobutil-add-output rampage_idr_bb "$rampage_idr_bb" --class=file
    dx-jobutil-add-output rampage_idr_png "$rampage_idr_png" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
