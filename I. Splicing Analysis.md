


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
- [VAST-TOOLS](https://github.com/vastgroup/vast-tools): Ben Blancoe's lab. Used by Raphaelle.
- [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html): focused on differential exon usage. [Vignette](http://127.0.0.1:12657/library/DEXSeq/doc/DEXSeq.pdf).
- JunctionSeq is like DEXSeq with junction reads included (and is written by the QoRTs team). JunctionSeq vignette - they have a great walkthrough that ... walks you through the whole process from beginning to end inc. QoRTs
- rMATS: useful for comparing with other ENCODE datasets
- MAJIQ is also good but parsing the output is a bit annoying (but the default was the best looking one of the lot!)
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
devtools::install_github("kcha/psiplot")
install.packages("optparse")
#to deactivate environment
> source deactivate
```

PSI = percent spliced in

Set up R with modules loaded in conda environment
```bash
# move INCLUSION_LEVELS_FULL-Hsa6-hg19.tab to local terminal (in bioinformatics/)
ml R/3.5.1-foss-2016b-bare
R
library("ggplot2")
library("optparse")
library("MASS")
library("RColorBrewer")
library("reshape2")
library("grid")
library("devtools")
library("psiplot")
library("optparse")
#quit R in cluster
q()
```
run vastools compare in conda environment outside of R
```bash
#run on cluster
IN=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/vast_out
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools compare $IN/INCLUSION_LEVELS_FULL-Hsa6-hg19.tab -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --plot_PSI -sp Hsa --GO"
```
Can use VAST-TOOLS here to calculate differentially expressed genes: `compare_expr`
Output file is created in directory of input file. This reports the differentially spliced AS events between the 2 groups (based on difference in average inclusion levels - delta PSI)

## Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#differential-splicing-analysis
Test for differential alternative splicing between 2 groups of samples.

```bash
ml Anaconda2
ml R

#create conda environment
source activate rtest

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
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools diff -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --sampleNameA=VCP --sampleNameB=CTRL -o $OUT -r 0.99 -m 0.2 -d diff.splicing_r0.99_mv0.2 -c 8"
# use the -m argument to set threshold for differential splicing events (default is 0.1)
```
Output file = diff.splicing.pdf & diff.splicing.tab

```bash
order tab file by MV value:
more diff.splicing_mv0.2.tab | sort -k6 -r | awk '{ if ($6 >= 0.2) { print } }' | awk '{ if ($5 >= 0) { print } }'

# count number of genes with MV > X:
awk '{ if ($6 >= 0.2) { print } }' diff.splicing_mv0.2.tab | wc -l
# count number of genes with +ve intron retention (VCP vs CTRL):
awk '{ if ($6 >= 0.2) { print } }' diff.splicing_mv0.2.tab | awk '{ if ($5 >= 0) { print } }' | wc -l
```
# DEXSeq

https://bioconductor.org/packages/release/bioc/html/DEXSeq.html
http://127.0.0.1:12657/library/DEXSeq/doc/DEXSeq.pdf
https://genome.cshlp.org/content/22/10/2008

Measures differential exon usage (DEU) which indicates alternative splicing. DEU also measures alternative transcript start sites & polyadenylation sites (differential usage of exons at 5' and 3' boundary of transcripts). 

Calculates ratio = number of transcripts from the gene containing this exon / total number of transcripts from the gene.
Then compares this ratio between conditions assessing the strength of the fluctuations (dispersion). 

Steps:
1. Alignment to genome using splice aware aligner eg STAR
2. Counts using DEX Seq Python scripts (utilises HTSeq)


```bash
ml HTSeq
ml Python/2.7.15-GCCcore-7.3.0-bare

# use 2 Python scripts (dexseq_prepare_annotation.py & dexseq_count.py)
# load R environment with DEXSeq in
source activate rtest
R
library("DEXSeq")

# Create output folder
mkdir -p DEXSeq

# set GTF - Ensembl (gencode)
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
# set GFF output
OUT=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_prepare_annotation.py

#run dexseq_prepare_annotation.py script
python $SCRIPT $GTF $OUT

#run count script
# set GFF input made in script 1
GFF=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_count.py
#set BAM input file
SAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*.sam
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/DEXSeq/SRR5483788.txt

#run dexseq_count.py
python $SCRIPT $GFF $SAM $OUT



for SAMPLE in $SAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	srun -N 1 -c 8 --mem=40GB --wrap="python $SCRIPT $GFF $SAMPLE ${SRRID}.txt"
done


#to deactivate environment
source deactivate
```


# Coverage for introns of interest
To perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.



# Perform GO analysis of gene list of differentially spliced genes

### Data Preparation

```bash
# convert splicing_diff.tab > .csv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < input.tab > output.csv
```
Open CSV in Excel
Copy the significant Gene Symbols into DAVID (see Gene Enrichment Chapter)

## topGO

### Prepare data from VAST tools splicing output

```r
library(tidyverse)
setwd("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/vast_tools/vast_out") #set working directory where splicing files are stored on cluster
diff_splicing <- read_tsv("diff.splicing_mv0.2.tab") #import differential splicing table .tab file
diff_splicing <- data.table(diff_splicing)

diff_splicing$pval <- as.numeric(diff_splicing$`E[dPsi]`)
diff_splicing$MV   <- as.numeric(diff_splicing$`MV[dPsi]_at_0.95`)

genelistUp <- factor( as.integer( diff_splicing$MV > .2 & diff_splicing$pval > 0) )
names(genelistUp) <- rownames(diff_splicing)
genelistDown <- factor( diff_splicing$MV > 0.2 & diff_splicing$pval < 0 )
names(genelistDown) <- rownames(diff_splicing)
```

### Run the topGO commands & make bar plots for BP, CC & MF (up & downregulated)

```r
### GO analysis
library(topGO)
library(GOstats)
library(goseq)
library(org.Mm.eg.db)
library(ggplot2)

### Test UPREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata      <- new( "topGOdata", ontology = "BP", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_BP_table   <- GenTable( myGOdata, goTestResults )
enrich.ur.bp  <- transform(UR_BP_table,result1 = as.numeric(result1))
enrich.ur.bp$value <- -log10(enrich.ur.bp$result1)
# Plot GO p-values as as bar plot
dat.ur.bp        <- enrich.ur.bp$value # the -log10(P-value)
names(dat.ur.bp) <- enrich.ur.bp$Term #the description of your GO term
barplot(height = dat.ur.bp, space = 0.5, horiz=T,las=1,font.size = 10)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.ur.cc <- transform(UR_CC_table,result1 = as.numeric(result1))
enrich.ur.cc$value <- -log10(enrich.ur.cc$result1)
# Plot GO p-values as as bar plot
dat.ur.cc        <- enrich.ur.cc$value # the -log10(P-value)
names(dat.ur.cc) <- enrich.ur.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.ur.cc,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_MF_table <- GenTable( myGOdata, goTestResults )
enrich.ur.mf <- transform(UR_MF_table,result1 = as.numeric(result1))
enrich.ur.mf$value <- -log10(enrich.ur.mf$result1)
### Plot GO p-values as as bar plot
dat.ur.mf        <- enrich.ur.mf$value # the -log10(P-value)
names(dat.ur.mf) <- enrich.ur.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.ur.mf,horiz=T,las=1, font.size = 20)

### TEST DOWNREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_BP_table <- GenTable( myGOdata, goTestResults )
enrich.dr.bp <- transform(DR_BP_table,result1 = as.numeric(result1))
enrich.dr.bp$value <- -log10(enrich.dr.bp$result1)
### Plot GO p-values as as bar plot
dat.dr.bp        <- enrich.dr.bp$value # the -log10(P-value)
names(dat.dr.bp) <- enrich.dr.bp$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.bp,horiz=T,las=1, font.size = 20)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.dr.cc <- transform(DR_CC_table,result1 = as.numeric(result1))
enrich.dr.cc$value <- -log10(enrich.dr.cc$result1)
### Plot GO p-values as as bar plot
dat.dr.cc        <- enrich.dr.cc$value # the -log10(P-value)
names(dat.dr.cc) <- enrich.dr.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.cc,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_MF_table <-GenTable( myGOdata, goTestResults )
enrich.dr.mf <- transform(DR_MF_table,result1 = as.numeric(result1))
enrich.dr.mf$value <- -log10(enrich.dr.mf$result1)
### Plot GO p-values as as bar plot
dat.dr.mf        <- enrich.dr.mf$value # the -log10(P-value)
names(dat.dr.mf) <- enrich.dr.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.mf,horiz=T,las=1, font.size = 20)
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTY1NDk5NTk4OCwtNzg1NDUzMDUxLDE1NT
A3OTE4ODYsMTI4MTcxNzIyNywxNjUyODQ5NjczLC0xMzczMTky
MjM3LC02NDY5Njg0MTUsLTEyMzkxODQ1ODksLTUwNzkwMDg2My
wtNjcwNjY3NDY1LDEwOTc5NzYwNzQsLTc1MzY5MzMyMCwtMTUw
MzcwNTE5OSwtMTI2OTg3NzE2MSwtMTI2NTg4MTY1OCwtMzA5Nj
EzMTUzLDc2NTkyNzg1NSwtMTQ3MDM3MjgxNiwtODMyMzEyNjA2
LC02NzAxNTM0MzhdfQ==
-->