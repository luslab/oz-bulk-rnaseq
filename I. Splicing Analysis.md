


> # Alternative Expression Analysis
https://www.nature.com/articles/nmeth.1503.pdf

Infer structural infromation about the transcript 
Infer the strand by examining splice site spanning reads
Each transcript isoform has very few exons & exon-exon junctions that are unique to that isoform.

Can ignore structure of full length transcript & focus on individual sequence features.

![enter image description here](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=info:doi/10.1371/journal.pcbi.1004393.g006)
yellow gene = noncoding RNA gene.
brown & green genes = coding genes
Few exon-exon spanning genes.

# Tools
- [VAST-TOOLS](https://github.com/vastgroup/vast-tools)
- rMATS- useful for comparing with other ENCODE datasets
- MAJIQ is also good but parsing the output is a bit annoying (but the default was the best looking one of the lot!)
- JunctionSeq is like DEXSeq with junction reads included (and is written by the QoRTs team). JunctionSeq vignette - they have a great walkthrough that ... walks you through the whole process from beginning to end inc. QoRTs
- Whippet is new and lightweight, but you can't really see what it is up to or the reads it has aligned (edited)
- MISO

# Gene Isoform counting

Gene isoforms are mRNA produced from the same locus but with different protein coding sequences.  5 modes of alternative splicing are recognised:
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

#run vast-tools on each FASTQ file separately. Dont specify output as all files need to be in same subfolder > output auto goes into a folder called vast_out. Run from the vast-tools directory
for SAMPLE in $FASTQ
do
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE"
done
```

## Merging Output
https://github.com/vastgroup/vast-tools#merging-outputs

If no technical replicates then skip this.

Merge the Aligned output files for technical replicates when read coverage for independent replicates is not deep enough for a complete AS analysis.  Ideally have >150 million reads per sample for VAST-TOOLS AS analysis. 

## Combining results

Combines aligned files that are stored in the folder `to_combine` to form one final table called `INCLUSION_LEVELS_FULL-Hsa6-hg19.tabz`. This is the file that you send to differential splicing command. Can specify hg38. The output directory contains the sub-folders to combine..
```bash
#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/vast_out/

#  create the old legacy version INCLUSION_TABLE.tab single output then specify `--noANNOT`
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa --noANNOT"

#check output
head INCLUSION_LEVELS_FULL-Hsa6-hg19.tab

# run vast-tools combine using new v2.0.0 ANNOT tool 
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa"
# This produces 5 INCLUSION_TABLE files in the raw_incl folder. 
##ANNOT = identifies & profiles annotated exons
##COMBI = splice site based module
##EXSK = 
##MIC = 
##MULTI =
```

## Compare Groups
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples

As CAMP R module doesnt have psiplots R package need to create a conda environment to install packages:
```bash
ml Anaconda2
#create conda environment
conda create -n rtest r-essentials r-devtools
source activate rtest
# install the package normally by calling R
R
install command
load package

#to deactivate environment
> source deactivate
```

PSI = percent spliced in
Output file is created in directory of input file. 

```bash

# move INCLUSION_LEVELS_FULL-Hsa6-hg19.tab to local terminal (in bioinformatics/)
R
library("ggplot2")
library("optparse")
library("MASS")
library("RColorBrewer")
library("reshape2")
library("grid")
library("psiplot")
#quit R in cluster
q()

IN=/Users/ziffo/bioinformatics
vast-tools compare $IN/INCLUSION_LEVELS_FULL-Hsa6-hg19.tab -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --plot_PSI -sp Hsa --GO

#run on cluster
IN=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/vast_out
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools compare $IN/INCLUSION_LEVELS_FULL-Hsa6-hg19.tab -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --plot_PSI -sp Hsa --GO"

```
Can use VAST-TOOLS here to calculate differentially expressed genes: `compare_expr`

## Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#differential-splicing-analysis
Test for differential alternative splicing between 2 groups of samples.

```bash
ml R
R
library("ggplot2")
library("optparse")
library("MASS")
library("RColorBrewer")
library("reshape2")
library("grid")
library("psiplot")
#quit R in cluster
q()

OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/vast_out

sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools diff -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --sampleNameA=VCP --sampleNameB=CTRL -o $OUT -d diff.splicing -c 8"
```

## Plotting

```bash
vast-tools plot outputdir/significant_events.tab
```

# Coverage for introns
To perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.






<!--stackedit_data:
eyJoaXN0b3J5IjpbMTczMzA0MTU5MiwyMDE3MDExMzAyLC00Mz
kxNjE3NTgsMTM3NDg4MTYxNiwxMTcyNjQ0Mjg1LDEyMzkxMjE4
MCwtODc2NzA3NjQ4LC0xNDI4MzAwNzYwLDYwMzc2NTQ5NCwtNT
M4MTMyMjA5LC0xNTQxNDAzMzczLDE3MDQ2MDk1NTAsMTc1MjI3
MDQ2NiwtMTE2NDE2OTU4NSwtODc1OTUzMDgxLC0xMzQ5ODAzMj
Y5LC0xNjQ3MzkyNDg0LDExMzA2MDQ2MDYsLTE5Njg0NjQ5OTks
MTc1ODk1NTU0MV19
-->