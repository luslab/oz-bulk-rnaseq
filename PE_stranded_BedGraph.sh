#!/bin/bash

# Usage: SE_stranded_TDF.sh -p [OUTDIR FOR TEMP FILES] -b [PATH TO BAM FILE] -n [PATH TO OUTPUT FILE]

# Example:
# cd /farm/home/luisie01/Riccio/Exp_1/Run_gap_1/utr_id
# /farm/home/luisie01/Scripts/3-UTR-pipeline/rat_study/PE_stranded_Coverage.sh -p


while getopts p:b:n: option
do
        case "${option}"
        in
                p) OUTDIR=${OPTARG};;
                b) BAM=${OPTARG};;
                n) OUT=${OPTARG};;
        esac
done


cd $OUTDIR
echo pwd

# 1. Calculate coverage separately for each strand

# Watson-strand:
genomeCoverageBed -ibam $BAM -bg -split -strand + > $OUT.neg.bedgraph
# Crick-strand:
genomeCoverageBed -ibam $BAM -bg -split -strand - > $OUT.pos.bedgraph

echo "finished reading coverage"


# 2. Modify all numbers that are bigger than 9999
/bin/ls -d *.bedgraph | while read Line
do
    echo $Line
    temp=$Line.temp
    /bin/sed 's/[0-9]*.[0-9]*e+0[1-9]$/99999/' $Line > $temp && mv $temp $Line
    echo " "
    echo " "
done
echo "finished modify numbers"


# 3. Convert (+) into (-) for negative strand
sed 's/\([0-9]*$\)/-\1/' $OUT.neg.bedgraph > $OUT.neg.modif.bedgraph

echo "finished conversion negative branch"

# 5. Merge positive and negative strand
cat $OUT.neg.modif.bedgraph $OUT.pos.bedgraph > $OUT.bedgraph
echo "finished merging files"

# 7. Sort bedgraph
igvtools sort $OUT.bedgraph $OUT.sorted.bedgraph
igvtools sort $OUT.neg.bedgraph $OUT.neg.sorted.bedgraph
igvtools sort $OUT.pos.bedgraph $OUT.pos.sorted.bedgraph
echo "finished sorting files"

# 8. Binarise for IGV
igvtools toTDF $OUT.sorted.bedgraph $OUT.tdf hg19.chrom.sizes
igvtools toTDF $OUT.pos.sorted.bedgraph $OUT.pos.tdf hg19.chrom.sizes
igvtools toTDF $OUT.neg.sorted.bedgraph $OUT.neg.tdf hg19.chrom.sizes
echo "finished toTDF files"

# 8. Clean the directory
#/bin/rm -rf m1.pos.bedgraph m2.pos.bedgraph m1.neg.bedgraph m2.neg.bedgraph *.bam

echo "I am done with creating per nucleotide coverage file"
