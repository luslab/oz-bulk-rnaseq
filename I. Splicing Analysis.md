


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
- [rMATS](http://rnaseq-mats.sourceforge.net/): useful for comparing with other ENCODE datasets
- [MISO](http://genes.mit.edu/burgelab/miso/)
- MAJIQ is also good but parsing the output is a bit annoying (but the default was the best looking one of the lot!)
- Whippet is new and lightweight, but you can't really see what it is up to or the reads it has aligned (edited)
- [LeafCutter](https://www.nature.com/articles/s41588-017-0004-9) https://github.com/davidaknowles/leafcutter
- [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#overview-of-alternative-splicing-workflow)

# Gene Isoform counting

Gene isoforms are mRNA produced from the same locus but with different protein coding sequences.  5 modes of alternative splicing are recognised:
1.  Exon skipping or cassette exon.
2.  Mutually exclusive exons.
3.  Alternative donor (5') site.
4.  Alternative acceptor (3') site.
5.  Intron retention.

![Alternative Splicing](https://en.wikipedia.org/wiki/Protein_isoform#/media/File:Alternative_splicing.jpg)

# VAST-TOOLS

## Approach using iPSC differentiation

VAST-tools output = any splicing changes OVER TIME. i.e not just retained introns
Focus is how splicing patterns changed over time (Day 0; Day 7; Day 14 & Day 21) in the CTRL & VCP groups. Rather than directly compare splicing differences at each stage between CTRL & VCP, I compared the groups of genes at each stage which were exhibiting changes in splicing over time in control group versus VCP group. 
Most changes in IR in CTRL occurred at day 14 whilst in VCP the same events were occurring at day 7, premature IR/splicing in VCP mutants.

The `Splicing_VASTOOLS.sh` script is located in `/home/camp/ziffo/working/oliver/scripts/intron_retention`

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
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools
# Not necessary to create output folder as VAST-tools auto creates an output folder within your current working directory named "vast_out" and within that is "tmp" and "to_combine" folders

#set shortcuts
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D0_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D7_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D14_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D21_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D35_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D112_samples/SRR*_1.fastq

#run vast-tools on each FASTQ file separately at each time point. Dont specify output as all files need to be in same subfolder > output auto goes into a folder called vast_out. Run from the vast-tools directory
for SAMPLE in $FASTQ
do
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE -sp Hsa"
done
```

## Merging Output
https://github.com/vastgroup/vast-tools#merging-outputs

Merge the aligned output files for technical replicates when read coverage for independent replicates is not deep enough for a complete AS analysis.  Ideally have >150 million reads per sample for VAST-TOOLS AS analysis.  Raphaelle merged the 3 samples for each of VCP & CTRL at each time point.
If no technical replicates then skip this.

```bash
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools
#Prepare config_file from sample group txt file (needs 2 columns: fastq file name & group separated by a tab)
awk '{print $3"\t"$2}' /home/camp/ziffo/working/oliver/scripts/intron_retention/VASTOOLS_merge_groups.txt | tail -31 > /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/config_file

CONFILE=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/config_file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out

#Run it -- takes about 10 minutes on the 31 samples
vast-tools merge --groups ${CONFILE} --o $OUT --sp Hsa --move_to_PARTS
```

## Combining results
ml R
https://github.com/vastgroup/vast-tools#combining-results

About 2G of memory required; completed in about 10 min
Combine aligned files that are stored in the folders `to_combine`  to form one final table called `INCLUSION_LEVELS_FULL-Hsa6-hg19.tabz` (if using old legacy version) or 5 tables in v.2.0.0. This is the table that you send to differential splicing command. Can specify hg38. The output directory contains the sub-folders to combine. 
```bash
#move to to_combine directory with merged aligned files
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out/to_combine
#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out

#  create the old legacy version INCLUSION_TABLE.tab single output then specify `--noANNOT`
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa -v --noANNOT"

# run vast-tools combine using new v2.0.0 ANNOT tool - identifies annotated exon-exon reads
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa"
# This produces 5 INCLUSION_TABLE files in the raw_incl folder. 
##ANNOT = identifies & profiles annotated exons. This is the novel feature of v.2.0.0
##COMBI = splice site based pipeline
##EXSK = only alternative Exons are considered. Transcript based pipeline, single.
##MIC = microexon pipeline
##MULTI = transcript based pipeline, multiexon

#check output
## old legacy table
cd vast_out
head INCLUSION_LEVELS_FULL-Hsa6-hg19.tab
## new v2.0.0
cd vast_out/raw_incl
head INCLUSION_LEVELS_ANNOT-Hsa14-n.tab
head INCLUSION_LEVELS_COMBI-Hsa14-n.tab
head INCLUSION_LEVELS_EXSK-Hsa14-n.tab
head INCLUSION_LEVELS_MIC-Hsa14-n.tab
head INCLUSION_LEVELS_MULTI-Hsa14-n.tab
```
Format of the combine output is as follows: 
https://github.com/vastgroup/vast-tools/blob/master/README.md#combine-output-format

## Compare Groups & Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples

As CAMP R module doesn't have psiplots R package need to create a conda environment to install packages. Once this environment has been made it is saved for future and can be activated using `source activate`
```bash
ml Anaconda2
source activate rtest

# To create a conda environment from scratch
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

There are 2 approaches to comparing: 
1. Time effect (delta change in differential splicing between different time points of iPSC differentiation)
2.  Mutant effect (VCP vs CTRL at specified time points)

### Time Effect

```bash
mkdir /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/time_effect
INFILE=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out/INCLUSION_LEVELS_FULL-Hsa14-hg19.tab
```

Run the script in `/home/camp/ziffo/working/oliver/scripts/intron_retention/Splicing_VASTOOLS.sh`
Run the Compare step for each of the WT & VCP time-point comparisons

### Mutant Effect



```bash
#run on cluster
IN=/home/camp/ziffo/working/oliver/projects/airals/splicing/vast_tools/vast_out/INCLUSION_LEVELS_FULL-Hsa6-hg19.tab
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools compare $IN -a SRR5483788_1,SRR5483789_1,SRR5483790_1 -b SRR5483794_1,SRR5483795_1,SRR5483796_1 --min_dPSI 15 --min_range 5 --GO -sp Hsa"
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

## Coverage for introns of interest
To perform the more focussed analysis on the 167 retained introns identified using VAST-tools, use this script in R which  calculates a ratio of intron sequence coverage and surrounding exons. 

1. First check that the gene where the event is occurring is expressed. 
2. Then check the event exhibits a change over time of at least 10%. 
3. Finally inspect visually all selected IR events occurring at Day 7 NPC stage in VCP mutant and at Day 14 pMN stage in CTRL to end up with the list of 167 events to get a high confidence list.

Import the results obtained from VAST-tools (remember analysis in VCP & CTRL over time performed initially independently) 
The `import_VARTOOLS.R` script is located in `/home/camp/ziffo/working/oliver/scripts/intron_retention`

To then get a value of IR across diverse data-sets I then wrote the custom code that computed the ratio between coverage of the intron versus average coverage of the neighbouring exons. Then select the events of interest. The `get_relative_coverage_interactive.R` script is located in `/home/camp/ziffo/working/oliver/scripts/intron_retention`

# DEXSeq

https://bioconductor.org/packages/release/bioc/html/DEXSeq.html
https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf
https://bioconductor.org/packages/release/bioc/manuals/DEXSeq/man/DEXSeq.pdf
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
ml R

# use 2 Python scripts (dexseq_prepare_annotation.py & dexseq_count.py)
# load R environment with DEXSeq in
source activate rtest
R
library("DEXSeq")
#exit R
q()

# Create output folder
mkdir -p DEXSeq

### RUN ANNOTATION SCRIPT
# set GTF - Ensembl (gencode)
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
# set GFF output
OUT=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_prepare_annotation.py

#run dexseq_prepare_annotation.py script
python $SCRIPT $GTF $OUT


### RUN COUNT SCRIPT
# set GFF input made in annotation script
GFF=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_count.py
#set BAM input file
SAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR5483796.sam
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/DEXSeq/SRR5483796.txt

#run dexseq_count.py
python $SCRIPT $GFF $SAM $OUT

#to deactivate environment
source deactivate
```
Now switch to R to run the DEXSeq analysis.

```r
# make new Rproject & save it in the relevant CAMP Cluster splicing folder.
setwd("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/DEXSeq")
suppressPackageStartupMessages( library( "DEXSeq" ) )

### PREPARE DATA ###
# read in files
countFiles = list.files("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/DEXSeq/", pattern=".txt", full.names=TRUE) 
basename(countFiles)
flattenedFile = list.files("/Volumes/lab-luscomben/working/oliver/genomes/annotation/", pattern="gff", full.names=TRUE) 
basename(flattenedFile)

# create sample table: 1 row for each library. columns for file name & read counts, sample
sampleTable = data.frame(
  row.names = c( "SRR5483788", "SRR5483789", "SRR5483790", "SRR5483794", "SRR5483795", "SRR5483796" ), 
  condition = c("VCP", "VCP", "VCP", "CTRL", "CTRL", "CTRL" ) )
# check table
sampleTable
# create DEXSeqDataSet 
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
# see structure of dxd
sampleAnnotation( dxd )
colData(dxd)
# see first 5 rows from count data - 12 columns (first 6 = reads mapping to exons, last 6 = counts mapping to rest of exons from same gene)
head( counts(dxd), 5 )
#show details of exon bins annotation
head( rowRanges(dxd), 3 )

### DEXSeq ANALYSIS ###
#normalise for different sequencing depths between samples
dxd = estimateSizeFactors( dxd )
# Perform differential exon usage test in a single command (this by passes the steps below - takes 10 mins)
dxd = DEXSeq(dxd,
       fullModel=design(dxd),
       reducedModel = ~ sample + exon,
       BPPARAM=MulticoreParam(workers=1),
       fitExpToVar="condition")

# OPTIONAL: dispersion estimation - takes 30 mins!
dxd = estimateDispersions( dxd )
# plot the per-exon dispersion estimates vs mean normalised count
plotDispEsts( dxd )

# test for differential exon usage (DEU) between conditions
dxd = testForDEU( dxd )
# estimate relative exon usage fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

### RESULTS ###
# summarise the resuls
dxr1 = DEXSeqResults( dxd )
# description of each column in DEXSeqResults
mcols(dxr1)$description
# how many exons are significant (with false discovery rate 10% - can change this threshold to more stringent)
table ( dxr1$padj < 0.1 )
# how many genes are affected
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
# MA plot: log fold change vs avg normalised count per exon. Significant exons in red
plotMA( dxr1, cex=0.8 )
### Detailed overview of analysis results in HTML - see the geneIDs of differentially spliced genes
DEXSeqHTML( dxr1, FDR=0.1, color=c("#FF000080", "#0000FF80") )

### VISUALISATION ###
# visualise results for individual genes
plotDEXSeq( dxr1, "ENSG00000116560", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# visualise transcript models to see isoform regulation
plotDEXSeq( dxr1, "ENSG00000116560", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# visualise individual samples rather than model effect estimates.
plotDEXSeq( dxr1, "ENSG00000116560", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# remove overall changes from the plots
plotDEXSeq( dxr1, "ENSG00000116560", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
```

# rMATS
ml Python/2.7.15-GCCcore-7.3.0-bare
ml numpy
ml OpenBLAS
ml STAR
ml ScaLAPACK
ml GSL
ml GCC

[User Tutorial](http://rnaseq-mats.sourceforge.net/user_guide.htm)
http://rnaseq-mats.sourceforge.net/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280593/

MATS is a computational tool to detect differential alternative splicing events from RNA-Seq data. The statistical model of MATS calculates the P-value and false discovery rate that the difference in the isoform ratio of a gene between 2 conditions exceeds a given user-defined threshold. From the RNA-Seq data, MATS can automatically detect and analyze alternative splicing events corresponding to all major types of alternative splicing patterns. MATS handles replicate RNA-Seq data from both paired and unpaired study design.

![enter image description here](http://rnaseq-mats.sourceforge.net/splicing.jpg)





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
eyJoaXN0b3J5IjpbLTQ1MDk4NTMyMywtNDE2MjAwNTUwLC0xND
A2NzQxNTMsMTc3Mzk0MjI5NywzNjI2Mzk1MTQsMTM5NTkyOTky
MiwtMTExOTMyMDM3NSwtMTg0NTIxMjQ0NywtMjA5MjAwMzgzNC
w3MDc5MzM0NjAsMTM3MjkxNzM4OCwtMTI1NzQyMzQwMiw4NjEz
OTQyNjEsLTE5MDY2ODU3MiwtMTYyMzU1NzA4LDM5NzYyODM4MC
wtNDg1ODcwMTQ0LC02MjYyNTI5NjcsLTEzMTgwMzA2NDcsOTQ0
ODY0NjAwXX0=
-->