


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

#run vast-tools on each FASTQ file separately. Dont specify output as all files need to be in same subfolder. Run from the vast-tools directory
for SAMPLE in $FASTQ
do
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE"
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
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/

# run vast-tools combine
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa"
```

## Compare Groups
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples

PSI = percent spliced in
Output file is created in directory of input file. 

```bash
IN=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools

sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools compare $IN/INCLUSION_LEVELS_FULL-Hsa6-hg19.tab -a SRR5483788,SRR5483789,SRR5483790 -b SRR5483794,SRR5483795,SRR5483796 --plot_PSI -sp Hsa --GO"
```
Can use VAST-TOOLS here to calculate differentially expressed genes: `compare_expr`

## Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#differential-splicing-analysis
Test for differential alternative splicing between 2 groups of samples.

```bash
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/

vast-tools diff -a $OUT/SRR5483788,$OUT/SRR5483789,$OUT/SRR5483790 -b $OUT/SRR5483788,$OUT/SRR5483789,$OUT/SRR5483790 --sampleNameA=VCP --sampleNameB=CTRL -o $OUT -d diff.splicing -c 8
```

## Plotting

```bash
vast-tools plot outputdir/significant_events.tab
```

# Coverage for introns
To perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.






<!--stackedit_data:
eyJoaXN0b3J5IjpbMTE0MzM3NDkyNCwxMjAxOTA0Nzc2LC0xNz
QxODEwNzM3LDEzMjMxMzkzNzcsNzg0NzcxNTksMTI4OTczMTk2
NiwtNjM0MDE1NTkxLC0xMTMwMzk2NTI3LC0xNjk1NzE5NzY2LD
E3Mzg4NTU4MTIsMTk2MjkwNDk5MywyMDU4NzQxNzA3LDYyMjQ2
ODkxNCwyNDE5ODMzODYsLTE4MTE4MzI4MTEsLTE3MjkwNTExOT
IsLTE2ODg0NDYxMzQsLTEwNTY5NTEyNzYsNzMxOTgzNzQ2LDUz
NDMwNTY4NF19
-->