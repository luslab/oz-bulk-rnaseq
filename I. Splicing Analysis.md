


> # Alternative Expression Analysis
https://www.nature.com/articles/nmeth.1503.pdf

By infer structural infromation about the transcript 
Infer the strand by examining splice site spanning reads
Each transcript isoform has very few exons & exon-exon junctions that are unique to that isoform.

Can ignore structure of full length transcript & focus on individual sequence features.

![enter image description here](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=info:doi/10.1371/journal.pcbi.1004393.g006)
yellow gene = noncoding RNA gene.
brown & green genes = coding genes
Few exon-exon spanning genes.

# Tools
[VAST-TOOLS](https://github.com/vastgroup/vast-tools)
ENCODE use rMATS
MAJIQ is also good and interesting but parsing the output used to be a bit annoying (but the default was the best looking one of the lot!) (edited)
JunctionSeq is like DEXSeq with junction reads included (and is written by the QoRTs team)
Whippet is new and lightweight, but you can't really see what it is up to or the reads it has aligned (edited)
There is also MISO.
Most recently I've been using rMATS as I have been comparing stuff with other ENCODE datasets
But previously I used DEXSeq and JunctionSeq as it sort of follows on methodologically from DEseq
But I have tried MAJIQ and Whippet too... but they get more complex in terms of splice graphs
Whippet is in Julia, I think MAJIQ is python
JunctionSeq vignette - they have a great walkthrough that ... walks you through the whole process from beginning to end inc. QoRTs
And they improve on some of the flaws/limitations of DEXSeq
And you can effectively run a DEXSeq analysis alongside to see the difference

# Gene Isoform counting

Gene isoforms are mRNA produced from the same locus but with diferent protein codeing sequences:

5 modes of alternative splicing are recognized:

1.  Exon skipping or cassette exon.
2.  Mutually exclusive exons.
3.  Alternative donor (5') site.
4.  Alternative acceptor (3') site.
5.  Intron retention.

![Alternative Splicing](https://en.wikipedia.org/wiki/Protein_isoform#/media/File:Alternative_splicing.jpg)

# VAST-TOOLS

There are 4 steps:
1. Alignment (this is redone in VAST-TOOLS using bowtie - needs loading)
2. Combine outputs into 1 summary table
3. Differential Splicing analysis
4. Plot the output

## Alignment
ml R
ml Bowtie

https://github.com/vastgroup/vast-tools#alignment

Use untrimmed fastq files - use raw reads. Define reference genome species (Hsa = human). 
```bash
# Create output folder
mkdir -p vast_tools

#set FASTQ input file
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D7_samples/SRR54837*_1.fastq

#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/

#run vast-tools on each FASTQ file separately
for SAMPLE in $FASTQ
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE -o ${OUT}/${SRRID}"
done
```

## Merging Output
https://github.com/vastgroup/vast-tools#merging-outputs

If no technical replicates then skip this.

Merge the Aligned output files for technical replicates when read coverage for independent replicates is not deep enough for a complete AS analysis.
Ideally have >150 million reads per sample for VAST-TOOLS AS analysis. 

## Combing results

Combines aligned files to form one single summary table. This is the file that you send to differential splicing command. Can specify hg38. The output directory contains the sub-folders to combine..
```bash
#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/SRR5483788/
vast-tools combine -o $OUT -sp Hsa

OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/SRR54837*/
for SAMPLE in $OUT
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o ${OUT}/${SRRID}/ -sp Hsa"
done
```
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
#set Counts.txt output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/featureCounts/D7_samples/featureCounts/

#run featureCounts on each BAM file separately
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=24GB --wrap="featureCounts -a $GTF -g gene_name -o ${OUT}_${SRRID} $SAMPLE"
done

## Compare Groups
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples

PSI = percent spliced in
Output file is created in directory of input file. 

```bash
PATH=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools

vast-tools compare INCLUSION_TABLE.tab -a $PATH/SRR5483788,$PATH/SRR5483789,$PATH/SRR5483790 -b $PATH/SRR5483794,$PATH/SRR5483795,$PATH/SRR5483796 --plot_PSI -sp Hsa --GO
```
Can use VAST-TOOLS here to calculate differentially expressed genes: `compare_expr`

## Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#differential-splicing-analysis
Test for differential alternative splicing between 2 groups of samples.

```bash
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/

vast-tools diff -a $PATH/SRR5483788,$PATH/SRR5483789,$PATH/SRR5483790 -b $PATH/SRR5483788,$PATH/SRR5483789,$PATH/SRR5483790 --sampleNameA=VCP --sampleNameB=CTRL -o $OUT -d diff.splicing -c 8
```

## Plotting

```bash
vast-tools plot outputdir/significant_events.tab
```

# Coverage for introns
To perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.






<!--stackedit_data:
eyJoaXN0b3J5IjpbODk1MDk1NTE0LC0xMTMwMzk2NTI3LC0xNj
k1NzE5NzY2LDE3Mzg4NTU4MTIsMTk2MjkwNDk5MywyMDU4NzQx
NzA3LDYyMjQ2ODkxNCwyNDE5ODMzODYsLTE4MTE4MzI4MTEsLT
E3MjkwNTExOTIsLTE2ODg0NDYxMzQsLTEwNTY5NTEyNzYsNzMx
OTgzNzQ2LDUzNDMwNTY4NCwtMTA1MTMzOTkyMCwtMTE0NjE4Nz
E3LC01NDIzMDgzNjldfQ==
-->