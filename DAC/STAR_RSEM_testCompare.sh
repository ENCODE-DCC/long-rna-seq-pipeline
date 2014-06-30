### compare test for STAR_RSEM pipeline
### usage: STAR_RSEM_testCompare.sh <testDir1> <testDir2>
### output: compared file names, no differences should be reported

echo Log.final.out
diff <(awk 'NR>4{print}' $1/Log.final.out) <(awk 'NR>4{print}' $2/Log.final.out) | head

echo Aligned.sortedByCoord.out.bam
diff  <(samtools view $1/Aligned.sortedByCoord.out.bam) <(samtools view $2/Aligned.sortedByCoord.out.bam) | head

echo Quant.isoforms.results
diff  <(cut -f1-8 $1/Quant.isoforms.results) <(cut -f1-8 $2/Quant.isoforms.results) | head
echo Quant.genes.results
diff  <(cut -f1-7 $1/Quant.genes.results) <(cut -f1-7 $2/Quant.genes.results)| head

for ii in `cd $1; ls *bw`
do
    echo $ii
    diff $1/$ii $2/$ii | head
done
