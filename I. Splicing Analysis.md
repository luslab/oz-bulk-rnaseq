


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
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools

#run vast-tools on each FASTQ file separately
for SAMPLE in $FASTQ
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE -o ${OUT}/${SRRID}"
done
```

## Merging Output
https://github.com/vastgroup/vast-tools#merging-outputs

Merge the Align output files for technical replicates when read coverage for independent replicates is not deep enough for a complete AS analysis.
Ideally have >150 million reads per sample for VAST-TOOLS AS analysis. 

## Combing results

Combines aligned files to form one single summary table. This is the file that you send to differential splicing command. Can specify hg38. The output directory contains the sub-folders to combine..
```bash
#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/
vast-tools combine -o $OUT -sp Hsa
```

## Compare Groups
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples

PSI = percent spliced in
Output file is created in directory of input file. 

PATH=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools
```bash

vast-tools compare INCLUSION_TABLE.tab -a $PATH/SRR5483788,$PATH/SRR5483789, $PATH/SRR5483790 -b $PATH/SRR5483794,$PATH/SRR5483795, $PATH/SRR5483796 --min_dPSI 25 --min_range 5 -sp Hsa --GO
```
Can use VAST-TOOLS here to calcuate differentially expressed genes: `compare_expr`

## Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#differential-splicing-analysis




**Preparing an annotation:**
To **assess differential expression of exons**, create an annotation file where overlapping exons of different isoforms are split before running featureCounts. Use `dexseq_prepare_annotation.py` script of DEXSeq package or `QoRTs`.

@Raphaelle used: 
 **Splicing analysis**  
   
-   I first used  [VAST-tools](https://github.com/vastgroup/vast-tools) which performs alignment for you. So basically you submit your fastq files directly. Have a look at the GitHub vignette as it is rather complete however please do not hesitate to contact me if you want help with shell scripting as you will need to run this as a loop. Or Nobby will certainly be happy to help on CAMP (I am working from UCL cluster and have never logged onto CAMP).
-   Then to perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.

  **3' UTR isoforms analysis**
-   I have written an entire pipeline for this which I can explain and share with you scripts when needed. But basically I extract genome-wide coverage using bedtools, then extract regions of continuous coverage along genome, then intersect these with Ensembl annotated regions, extend 3' UTR. Finally to annotate all alternative 3' UTR isoforms I then run an algo which identifies shifts in coverage along 3' UTR which are expected to occur at PAS sites.

**SVD (singular value decomposition) analysis**

-   For doing this you can use the gene-level count table obtained from Kallisto. I wrote everything in R and I can send you some litterature which explains a bit the underlying math and idea. Also happy to speak about it over skype.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTUyOTc1NDY4OSw2MjI0Njg5MTQsMjQxOT
gzMzg2LC0xODExODMyODExLC0xNzI5MDUxMTkyLC0xNjg4NDQ2
MTM0LC0xMDU2OTUxMjc2LDczMTk4Mzc0Niw1MzQzMDU2ODQsLT
EwNTEzMzk5MjAsLTExNDYxODcxNywtNTQyMzA4MzY5XX0=
-->