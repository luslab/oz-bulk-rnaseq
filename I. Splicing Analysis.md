


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

# Coverage for introns
To perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.



# Perform GO analysis of gene list of differentially spliced genes



```r
### GO analysis
library(topGO)
library(GOstats)
library(goseq)
library(org.Mm.eg.db)
# subset results table to only genes with sufficient read coverage
resTested <- resLFC1[ !is.na(resLFC1$padj), ]
# remove decimal string & everything that follows from row.names (ENSEMBL) and call it "tmp": https://www.biostars.org/p/178726/ (alternatively use biomaRt )
temp=gsub("\\..*","",row.names(resTested))
#set the temp as the new rownames
rownames(resTested) <- temp
# construct binary variable for UP and DOWN regulated
genelistUp <- factor( as.integer( resTested$padj < .1 & resTested$log2FoldChange > 0 ) )
names(genelistUp) <- rownames(resTested)
genelistDown <- factor( as.integer( resTested$padj < .1 & resTested$log2FoldChange < 0 ) )
names(genelistDown) <- rownames(resTested)

### Test UPREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )
#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )
#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )

### TEST DOWNREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_BP_table <- GenTable( myGOdata, goTestResults )
#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )
#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
GenTable( myGOdata, goTestResults )

### Plot GO p-values as as bar plot
library("GOplot")

barplot(UR_BP_table, drop=TRUE, font.size = 10,showCategory=10000, title="Up Regulated Biological Process")


UR_BP_table$result1
table(UR_BP_table$Term)

UR_BP_table$value <- (UR_BP_table$result1)+1

pvalue <- UR_BP_table$result1
na.omit(pvalue)
log10(pvalue)

barplot(table(UR_BP_table$result1), names.arg=Term, main="Up Regulated Biological Process")  


height <- log10(result1)


circ <- circle_dat(UR_BP_table, sorted.df)
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTQ5ODIyNTM0OSwxNDEyNjk1MDE3LDIxND
UwNzY4MTgsLTMzNjg0MjI2OCwyMDAzNTA2NDAxLC0xNTMzNjI1
OTA0LC0xMTA1NjgyNzgsLTE3NDEwMjczODksNDI1MDU0NTQ4LC
0xNzQxMDI3Mzg5LDg1MDMxMDAwMCwtMTE2MjA3Mzk1LC0xNTgz
OTk0OTY0LC00Mzk5MDg1NTQsLTIxMzM2OTIwNDYsMTU3MTI3Mz
M2NywzOTYzMjExMjksMTczMzA0MTU5MiwyMDE3MDExMzAyLC00
MzkxNjE3NThdfQ==
-->