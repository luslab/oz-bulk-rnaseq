#!/bin/bash

# Usage: /farm/home/luisie01/Scripts/3-UTR-pipeline/rat_study/PE_stranded_Coverage.sh -p [OUTDIR FOR TEMP FILES] -b [PATH TO BAM FILE] -f [PROTOCOL USED: 1 IF ++/-- AND 2 OTHERWISE] -n [PATH TO OUTPUT FILE] -c [CHROMOSOME SIZE]



while getopts p:b:f:n:c: option
do
        case "${option}"
        in
                p) OUTDIR=${OPTARG};;
                b) BAM=${OPTARG};;
                f) PROTOCOL=${OPTARG};;
                n) OUT=${OPTARG};;
        esac
done


#CHROM=/farm/home/luisie01/Conservative_reference_annotation_rat/rn5/genome/rn5.chrom.sizes

# 1. Separate mates
if [ -f $OUTDIR/m1.bam ]; then
    echo No need to create BAM files
else
    echo "Create BAM files"
    /farm/babs/redhat6/bin/samtools view -b -f 0x0040 -h -o $OUTDIR/m1.bam $BAM
    /farm/babs/redhat6/bin/samtools view -b -f 0x0080 -h -o $OUTDIR/m2.bam $BAM
fi
echo "finished separation between left and right reads"

cd $OUTDIR

# 2. Calculate coverage separately for each strand

###
# M1
###

# Watson-strand:
/farm/home/luisie01/bin/genomeCoverageBed -bg -split -ibam m1.bam -strand + > m1.pos.bedgraph
# Crick-strand:
/farm/home/luisie01/bin/genomeCoverageBed -bg -split -ibam m1.bam  -strand - > m1.neg.bedgraph
echo "finished reading coverage m1"

###
# M2
###

# Watson-strand:
/farm/home/luisie01/bin/genomeCoverageBed -bg -split -ibam m2.bam  -strand + > m2.pos.bedgraph
# Crick-strand:
/farm/home/luisie01/bin/genomeCoverageBed -bg -split -ibam m2.bam  -strand - > m2.neg.bedgraph
echo "finished reading coverage m2"


# 3. Modify all numbers that are bigger than 9999
/bin/ls -d *.bedgraph | while read Line
do
    echo $Line
    temp=$Line.temp
    /bin/sed 's/[0-9]*.[0-9]*e+0[1-9]$/99999/' $Line > $temp && mv $temp $Line
    echo " "
    echo " "
done
echo "finished modify numbers"

# 4. Merge m1.pos with m2.neg and m1.neg with m2.pos
#    PROTOCOL=1 --> ++,--; PROTOCOL=2 --> +-,-+
if [ $PROTOCOL != 2 ]; then
    /bin/cat m1.pos.bedgraph m2.neg.bedgraph > $OUT.pos.bedgraph
    /bin/cat m1.neg.bedgraph m2.pos.bedgraph > $OUT.neg.bedgraph
else
    /bin/cat m1.neg.bedgraph m2.pos.bedgraph > $OUT.pos.bedgraph
    /bin/cat m1.pos.bedgraph m2.neg.bedgraph > $OUT.neg.bedgraph
fi

echo "finished merging files"

# 5. Convert (+) into (-) for negative strand
/bin/sed 's/\([0-9]*$\)/-\1/' $OUT.neg.bedgraph > $OUT.neg.modif.bedgraph

echo "finished conversion negative branch"

# 6. Merge positive and negative strand
/bin/cat $OUT.neg.modif.bedgraph $OUT.pos.bedgraph > $OUT.bedgraph
echo "finished merging files"

# 7. Sort bedgraph
/farm/home/luisie01/bin/IGVTools/igvtools sort $OUT.bedgraph $OUT.sorted.bedgraph
/farm/home/luisie01/bin/IGVTools/igvtools sort $OUT.neg.bedgraph $OUT.neg.sorted.bedgraph
/farm/home/luisie01/bin/IGVTools/igvtools sort $OUT.pos.bedgraph $OUT.pos.sorted.bedgraph
echo "finished sorting files"

# 8. Binarise for IGV
/farm/home/luisie01/bin/IGVTools/igvtools toTDF $OUT.sorted.bedgraph $OUT.tdf /farm/home/luisie01/bin/IGVTools/genomes/rn5.chrom.sizes
/farm/home/luisie01/bin/IGVTools/igvtools toTDF $OUT.pos.sorted.bedgraph $OUT.pos.tdf /farm/home/luisie01/bin/IGVTools/genomes/rn5.chrom.sizes
/farm/home/luisie01/bin/IGVTools/igvtools toTDF $OUT.neg.sorted.bedgraph $OUT.neg.tdf /farm/home/luisie01/bin/IGVTools/genomes/rn5.chrom.sizes
echo "finished toTDF files"

# 8. Clean the directory
#/bin/rm -rf m1.pos.bedgraph m2.pos.bedgraph m1.neg.bedgraph m2.neg.bedgraph *.bam

echo "I am done with creating per nucleotide coverage file"
