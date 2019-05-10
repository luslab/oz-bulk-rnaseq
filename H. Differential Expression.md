> # Differential Gene Expression Analysis

This covers downstream interpretation of expression & differential estimates. To compare two conditions look at the fraction of transcripts assigned to a specific gene over the total number of reads (total read number differs drastically between samples). The number of sequenced reads depends on:
1. expression level
2. read length
3. sequencing depth
4. expression of all other genes in the sample

# Normalisation

There are gene counting tools that normalise for sequencing depth at Read Quantification step, however DE analysis (DESeq & edgeR) normalise whilst doing DGE analysis. 

**Normalisation** & **Log Transforming** Read Counts is performed to ensure that systematic effects not related to the biological differences between samples are removed. 

2 reasons we need to normalise:
1. Differences in **Library sizes (sequencing depth)** between samples - RPKM, FPKM, TPM, CPM addresses this.
![enter image description here](https://lh3.googleusercontent.com/9qxYyeRJPOD_Iv8MRfRfbBbHVkUvpW-H9JQTCRYhdF_uYqRbTFI8ELV7OQfSeOcz_E6ea5bObIi9FA)
2. Differences in **Library composition (different sample type** e.g. neurones vs glia )
![enter image description here](https://lh3.googleusercontent.com/RcxEg6_OW73xGL4z_PDorRNJV9WcsHJoOmB09ZhHd4YTFfiQ6bgzXoRL2QsLgssjblzhXI9FIlJWWg)

## Normalisation Methods

- FPKM (RPKM): Reads/Fragments per kilobase of transcript per millions of read mapped. Fragment refers to read pairs from paired-end reads (counting fragments and not individual reads. 

FPKM normalises for gene size & sequencing depth using: **FPKM = (10^9^ * C ) / (N * L)**
C = number of mappable reads (fragments) for gene/transcript/exon etc
N = total number of mappable reads in the library
L = number of base pairs in the gene/transcript/exon etc (i.e. the size of gene length) 
FLAWED STATISTICALLY. These has been superseded by:
- Transcripts per Million (TPM); Counts/million (CPM): also normalises sequencing depth & gene length but in the reverse order. [TPM has superseded FPKM as it is a fairer statistical measure](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/).  Per Million helps keep the value sensible (avoids huge number of decimals) 
![enter image description here](https://lh3.googleusercontent.com/wB60H7Ih8MZmQk_7TAkj_Kc9KKO-vKpyvszwsf2xmpwOdV7q-QXqyGHkO_jDbczcHW2id7g1VW9w6w)
- TMM (trimmed mean of M values)
- upper quartile

# Filtering low read counts
[StatQuest video](https://www.youtube.com/watch?v=Gi0JdrxRq5s&feature=youtu.be) 
aka Independent Filtering or Mitigating the multiple testing problem
False positives are a huge problem when significance checking every gene in the genome: 5% of 20,000 statistical tests = 1000 false positives. 
FDR & Benjamini Hochberg method partly compensate for this BUT this adjustment looses some true positives > false negatives. With each statistical test we reduce the number of true positive p-values that survive the FDR adjustment.
![enter image description here](https://lh3.googleusercontent.com/y9gZel8cdmsYicJXqpgh7dP5erLxMosgCbl4C1_P3Z4jpYzhvwDplbTsi3q4hn1_2PVxtcOljdH3Iw)
To filter bogus tests, edgeR & DESeq2 use alternative filters:
- Very low read counts are not informative so they are removed > reduce number of tests
	- EdgeR removes all genes except thos with >1CPM in 2 or more samples. Note that CPM scaling factor is susceptible to sequencing depth when it is very high (loose relevant but lower read counts) or very low (e.g. sincle cell RNASeq, includes irrelevant low read counts)
![enter image description here](https://lh3.googleusercontent.com/LnuncxJbacK6DUrC12OPFhEFMeAS_1m9mpzzOfA2gme_6BL4ShdfnmYy_a5vQPObuR5B4yyzDoWETg)

Changing the CPM threshold to exclude genes drastically changes the number of signficant genes.
Red line = CPM 1 threshold: to strict - removes important genes
Blue line = CPM 0.2 threshold is better as it includes all the approapriate genes based on the curve (from peak)
![enter image description here](https://lh3.googleusercontent.com/caQjY78DDqLrBeG3DQDaV6wWwFnUVVWTphhmEZUP_QFl9bnKbdyvauR8gerHZjvQiy_qhYRH2yPtww)
DESeq2 calculates p-values before trying different CPM thresholds

## EdgeR vs DESeq2 Filtering
1. EdgeR looks at individual samples. Removes genes with <2 samples of CPM 1 or more.
2. DESeq2 looks at average normalised reads across all samples. If average is above CPM threshold then it keeps the gene.  This would be susceptible to huge outliers e.g.
![enter image description here](https://lh3.googleusercontent.com/2gKnrHbFKwU-7wMItYGJIf-NiI-h0JAlQ3o0TykWS8bXHrMdmxAkvHCiiEjXz8bbi8vp8YrDXHvlQQ)
To get around this DESeq2 has an outlier detection method used when there is >2 samples.
Note peak is similar with simi
![enter image description here](https://lh3.googleusercontent.com/7rmKBF-BM_NC1D80xxn0PrS-x8ggFZc7xwQEfJEHTz1qeZUxjhpi_fSBZfMaheDsbEOjPlbSp8LEdA)

### SVD_analysis.Rmd
Raphaelle normalises using the **SVD_analysis.Rmd** markdown. Uses a filtering function that is a 2 component mixture model to separate low vs high counts. Run `/Volumes/lab-luscomben/working/oliver/scripts/intron_retention/SVD analysis.Rmd` script directly in R studio

1. Import HTSeq raw gene count matrices
2. Gene annotation from ensembl
3. Label samples
4. Create colour palatte
5. Filter low count reads out
6. Normalise between samples using Limma

## Comparison Types

1.  _Within-sample_  comparisons: compare the expression of genes within the same experiment eg in this experiment does gene A express at a higher level than gene B?
2.  _Between-sample_  comparisons: compare the expression of a gene between experimental conditions aka pairwise comparison e.g. has the gene expression for gene A  changed across different experimental conditions? if yes, then differentially expressed (DE).

## Goals of DE analysis
The 2 tasks of DGE are to:
1. Estimate the magnitude **fold change** of differential expression in read counts between different conditions, accounting for differences in sequencing depth & variability
2. Estimate the **significance** of the difference, accounting for multiple testing

Null hypothesis = the mean read counts of genes in different samples are equal between different conditions.

Need to estimate the **mean** and **dispersion** from the read counts:
- precision depends on the number & variation of replicates (i.e. the statistical power). 
- For RNA-seq typically there is only 2-3 replicates creating poor precision. There are tools to compensate for this by looking across genes with similar expression to reduce a given genes variance towards the regressed values.

## Tools for DGE analysis

 - DESeq2 - Nobby uses this.
 - edgeR - better for false positives, less conservative, recommended if <12 replicates
 - ballgown
 - limma: raphaelles Rscript uses this (normalizeQuantiles function)
 - Cuffdiff is slow, cant support multifactored experiments, can detect differential isoforms, high false positives

![enter image description here](https://lh3.googleusercontent.com/LVvCl3GXhNzUx5lyTrHsr0z_ZmI0nb51TBiY1-53VifMuYW8HR9-X54sfLwoH5gFyqahHOm8_QaWhg "Comparison of DGE programs")

DESeq2 and edgeR are both good. Stick to one. https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/
See the [StatQuest videos](https://statquest.org/video-index/) on high-throughput sequencing analysis 

All these tools require a similar input = a matrix of counts: columns represent different samples; rows represent different genes; and the integers populating the matrix represent the counts of reads (if single end) or fragments (if paired end).

**Below are walkthroughs of this analysis using DESeq2, EdgeR and Ballgown.**

DESEq2 nor EdgeR use RPKM, FPKM, TPM etc. Both DESeq2 & EdgeR create a scaling factor for each sample that accounts for read depth & library composition.

# DESeq2
[StatQuest video on DESeq2 normalisation](https://www.youtube.com/watch?v=UFB993xufUU)
1. **Log^e^** of all gene count values. 0 counts become -Infinity. This allows normalisation between different tissue types & smooths effect of outlier read counts
2. **Average each Gene (row) logs** (i.e. the Geometric Average). This removes impact of outliers. Averaging a row with a -Inf value becomes -Inf. 
3. **Filter out** Genes with Infinity i.e. if any of the samples had a 0 count.
4. Subtract the average log value from log(counts). This is the same as log (Sample Count / Average Count) i.e. the **ratio of reads to average** of all samples
![enter image description here](https://lh3.googleusercontent.com/YcuQUEulQ5OEhaSOsx9pJwwuusIxW8DAJziJzdYTbg9aKA10309WTgo6-ca_q0GtdgR3AXS7gObE6g)
5. Calculate the **median of the ratios** for each sample. This further removes influence of outliers.
![enter image description here](https://lh3.googleusercontent.com/zmpaRfyzX8-2Ake2YMQpnpSkXV-vr87uNvIHfKlcp40Y2ajJI0cNsMYfGbSL-0IcNaqL5yE-SefMLg)
6. Convert medians to "normal numbers" to get the **Scaling Factors** for each sample. e ^sample_median^  
7. Divide the original read counts by the Scaling Factors
![enter image description here](https://lh3.googleusercontent.com/64nhKIz6HkNRl5TzsK_aafA-OocPa8DQiRJIJoVsWgnvNfrmhFR4dy-c3EBGg4qsYQPsfyRQHlQOCw)

All DESeq2 information is available in the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis)

DESeq takes counts --> read them into R --> normalise for sequencing depth differences
```r
#######################
# Load Data into R #
#######################

#Set working directory in R where output will go:
working_dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/counts"
setwd(working_dir)
#can also set the working directory in R manually: Session > Set Working Directory > Choose Folder in the Cluster to set as wd.

#Read in gene mapping
mapping=read.table("/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/htseq/ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1)
# Read in count matrix
rawdata=read.table("/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/htseq/htseq_counts_table.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
# Check dimensions
dim(rawdata)
# Require at least 25% of samples to have count > 25
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)
```

```r
# DESeq( ) function wraps around the following 3 functions:
DESeq.ds = estimateSizeFactors(DESeq.ds) #sequencing depth normalisation
DESeq.ds = estimateDispersions(DESeq.ds) #gene wise dispersion estimates
DESeq.ds = nbinomWaldTest(DESeq.ds) #fits a negative binomial GLM & applies Wald stats to each gene

DGE. results = results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05) #results function allows extraction of base means across sample, SEs, and basic stats for each gene.
summary(DGE.results)

# The DESeq Result behaves as a data frame:
head (DGE. results)
table (DGE. results $ padj < 0.05)
rownames ( subset (DGE. results , padj < 0.05) )

DSeqDataSet #contains all the information in R. 
rowRanges( ) # rows = variables (genes, transcripts, exons)  - has info about genes (chr, start, end, strand, ID)
rownames #is the unique sample names
colData( ) #columns = samples 
assay( ) #stores count data with genes and samples. similar to `countData`
```

## Using QoRTs Output

Perform DE analysis with DESeq2 in R

Write the Decoder text files: decoder.by.UID.txt & decoder.bySample.txt (identical to that used in QC of Aligned Reads chapter):
```
sample.ID	group.ID  qc.data.dir
SRR5483788	VCP  	QoRTs_SRR5483788
SRR5483789	VCP  	QoRTs_SRR5483789
SRR5483790	VCP  	QoRTs_SRR5483790
SRR5483794	CTRL 	QoRTs_SRR5483794
SRR5483795	CTRL 	QoRTs_SRR5483795
SRR5483796	CTRL 	QoRTs_SRR5483796
```
If there are technical replicates then merge them at this point. QoRTs allows count data to be combined across technical replicates. See step 4 (chapter 9, page 15) http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf 
If there are no technical replicates (as with Nat Comms paper) then skip this step.

```bash
#set QoRTS QC input
QC=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/QoRTs_SRR5483788/,./QoRTs_SRR5483789/,./QoRTs_SRR5483790/,./QoRTs_SRR5483794/,./QoRTs_SRR5483795/,./QoRTs_SRR5483796/
#set output directory
OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/QoRTs_counts

#run QoRTs command for each QC file
java -jar $EBROOTQORTS/QoRTs.jar mergeCounts --mergeFiles $QC --verbose $OUT
```
Run QoRTS>DESeq2 analysis in R
```r
#Set working directory where counts expression files exist
working_dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2"
setwd(working_dir)

# Load libraries
library(QoRTs)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library("vsn")
library("pheatmap")
library("ggrepel")
library("apeglm")
library("AnnotationDbi")
library("Homo.sapiens")
library("Glimma")
library("ReportingTools")
library("regionReport")
library(topGO)
library(GOstats)
library(org.Mm.eg.db)

res <- read.qc.results.data("/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/", 
                            decoder.files = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.byUID.txt", 
                            calc.DESeq2 = TRUE, calc.edgeR = TRUE); 

### Build the QC plots.
# The makeMultiPlot.all can be used to automatically generate a full battery of multi-plot figures (as png): 
makeMultiPlot.all(res, 
                  outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPlots/", 
                  plot.device.name = "png"); 
# As PDF: QoRTs offers multi-page pdf reports as an alternative, simply by using the plot.device.name parameter: 
makeMultiPlot.all(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPDFs/",
                  plot.device.name = "pdf"); 
# Print all the basic plots as seperate pngs: 
makeMultiPlot.basic(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/basicPlots/", 
                    separatePlots = TRUE);

# Extract size factors. QoRTs generates these to normalise all samples to a comparable scale allowing downstream comparison with DESeq2 or edgeR
get.size.factors(res, outfile = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/sizeFactors.GEO.txt");

### Run DESeq2 DE analysis
decoder.bySample <- read.table("/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.bySample.txt", 
                               header=T,stringsAsFactors=F); 
directory <- "/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC"; 
sampleFiles <- paste0(decoder.bySample$qc.data.dir, "/QC.geneCounts.formatted.for.DESeq.txt.gz" ); 
sampleCondition <- decoder.bySample$group.ID; 
sampleName <- decoder.bySample$sample.ID; 
sampleTable <- data.frame(sampleName = sampleName, fileName = sampleFiles, condition = sampleCondition); 
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition); 

# Run & Print DESeq2 results
dds
dds <- DESeq(dds);
res <- results(dds);
res
summary(res)
resLFC <- lfcShrink(dds, coef="condition_VCP_vs_CTRL", type="apeglm")
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)



# Sort the results data frame by the padj and foldChange columns. 
sorted = res[with(res, order(padj,  -log2FoldChange)),  ]  
# Turn it into a dataframe to have proper column names. 
sorted.df = data.frame("id"=rownames(sorted),sorted)
#write the table out:
write.table(sorted.df, file = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2/DESeq2.results.txt", sep="\t", col.names=NA, quote=FALSE);

### Normalise Counts
# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)
# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)
# Save the normalize data matrix.
write.table(dt, file="/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2/norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

### Data Transformation
# extract transformed values using vst (rapid) or rlog (slower)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
# extract the matrix of normalised values
head(assay(vsd), 6)
colData(vsd)
head(assay(rld), 6)
# Visually assess effect of transformation on variance. Plot SD of transformed data vs mean
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
```
The DESeq2.results.txt file describes changes between the 2 conditions e.g.

```bash
id             baseMean   baseMeanA     baseMeanB   foldChange  log2FoldChange    pval       padj
ERCC-00130      29681        10455        48907        4.67        2.22         1.16e-88    9.10e-87
ERCC-00108        808          264         1352        5.10        2.35         2.40e-62    9.39e-61
ERCC-00136       1898          615         3180        5.16        2.36         2.80e-58    7.30e-57
```
-   `id`: Gene or transcript name that the differential expression is computed for
-   `baseMean`: The average normalized value across all samples,
-   `baseMeanA`,  `baseMeanB`: The average normalized gene expression for each condition,
-   `foldChange`: The ratio  `baseMeanB/baseMeanA`,
-   `log2FoldChange`: log2 transform of  `foldChange`. When we apply a 2-based logarithm the values become symmetrical around 0. A log2 fold change of 1 means a doubling of the expression level, a log2 fold change of -1 shows show a halving of the expression level.
-   `pval`: The probability that this effect is observed by chance. Only use this value if you selected the target gene a priori.
-   `padj`: The adjusted probability that this effect is observed by chance. Adjusted for multiple testing errors.

```bash
#Sort by gene ID select only columns foldchange and log2FoldChange. The DESeq2.results.txt file is already sorted according to padj
cat DESeq2.results.txt | sort | cut -f 1,5,6 > DESeq2.table

#How many genes are significantly differentially expressed (i.e. padj < 0.05 in column 8)?
cat DESeq2.results.txt | awk ' $8 < 0.05 { print $0 }' > DESeq2.diffgenes.txt

#How many differentially expressed genes do we have?
cat DESeq2.diffgenes.txt | wc -l
```

## Using featureCounts Output
```bash
ml R

# Using featureCounts output (counts.txt) print out only columns representing Gene ID & sample abundances (ie remove intermediate columns - Chr, Start, End, Strand & length)
cat counts.txt | cut -f 1,7-14 > sample_counts.txt

#pass simple_counts.txt through the script specifying the design of the experiment 
## in this case = 3 x 3 (3 cases, 3 controls)
cat sample_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.txt
```
```R
# Load the library.
library(DESeq2)
library(magrittr)

# cat counts.txt | Rscript deseq2.r
# Produces a table with differentially expressed genes. on the standard output.
# To install the requirements run the program with the 'install` parameter.
# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Experimental design must be specified as: NxM at the command line", call.=FALSE)
}
first = args[1]
if (first == 'install') {
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    stop("Installation completed", call.=FALSE)}
# Extract the experimental design from the command line.
design = unlist(strsplit(first, 'x'))
# Find the desing counts.
cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])
# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)
# Read the data from the standard input.
countData = read.table("stdin", header=TRUE, sep="\t", row.names=1 )
# Build the dataframe from the conditions
samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)
#coldata = read.table(coldata_file, header=TRUE, sep="\t", row.names=1 )
# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
#Set the reference to be compared
dds$condition = relevel(dds$condition,"cond1")
# Run deseq
dds = DESeq(dds)
# Format the results.
res = results(dds)
# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]
# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)
# Write the table out.
write.table(sorted.df, file="", sep="\t", col.names=NA, quote=FALSE)
# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)
# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)
# Save the normalize data matrix.
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)
```

### Using Kallisto Output
```bash
ml R

# Using kallisto output (counts.txt) first remove empty rows (represent transcripts that were not detected at all) by keeping rows only where the sum of all reads is > 25
cat counts.txt | awk ' ($2+$3+$4+$5+$6+$7+$8+$9) > 25 { print $0 } ' > temp.txt

# Convert effective lengths to integers using ask & %3.0f formating operator (DESeq2 only accepts integers) - this is arranged for 3x3.
cat temp.txt | awk ' { printf("%s\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\n", $1, $2, $3, $4, $5, $6, $7) } ' > valid.txt
## for 4x4 samples
cat temp.txt | awk ' { printf("%s\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9) } ' > valid.txt

# Run the differential expression estimator
cat valid.txt | Rscript deseq1.r 3x3 > D7-results.txt
```

# edgeR
[StatQuest video on Edge R normalisation](https://www.youtube.com/watch?v=Wdt6jdi-NQo)
edgeR is a bioconductor package designed for DE of raw counts (it does not use RPKM, TPM etc).
1. Remove all untranscribed genes (0 counts)
2. Choose the reference sample to normalise all others against. A good reference sample should be the most average sample 
	2a. Scale each sample by total read counts
	2b. For each sample determine the value such that 75% of the scaled data are <= that it.
	2c. Average the 75th quantiles between samples
	2d. choose the sample whose 75th quintile is closest to the average
![enter image description here](https://lh3.googleusercontent.com/6R3NhK9CMdTjyNGj85TJCcQldoW1y0SE2O6hTUrqr0o0oV98t7y_2e4a-3On4lhJZpDnxV3m758DoQ)
3. Create a **table of log ratios to identify biased genes**. Do this separately for each sample relative to the reference. 
		3a. filter out biased genes (very up/down regulated genes) by using log fold differences
		3b. Compare remaining genes with Reference gene counts. Log2(Reference Gene/Sample Gene). If Reference gene count or sample Gene count = 0 then value is -Inf.
		3c. Remove gene with Inf values (i.e. no reads mapped to gene in 1 or both samples) 
![enter image description here](https://lh3.googleusercontent.com/a8vrZVwNslMuV6rGYX46XhqRFNB-ei7GC2QL4Qt6TIFm5bxYKN01sqZ9yV43zOt48Wmv_WfCZvjmDQ)
4. Calculate a **table to identify highly/lowly transcribed genes in across samples (Mean of Logs)**
	4a. Calculate Geometric mean for each gene Log. 
	4b. Remove Inf values again (no reads mapped to them in 1 or more samples).
	
5. Sort both Table of Log Ratios & Mean of Logs tables by Low to High.
	5a. Filter out top 30% & bottom 30% biased genes (Log Ratios)
	5b. Filter out top 5% & bottom 5% of the highly/lowly transcribed genes (Mean of Logs)

6. Use remaining genes to calculate **Scaling Factor**
	6a. Calculate weight average of remaining Log2 ratios (Weighted trimmed mean of log2 ratios). Genes with more reads have more weight.
	6b. Convert weighted average of of Log2 ratios to "normal numbers"
	6c. Centre the scaling factors around 1 (divide each raw scaling factor by their geometric mean)
![enter image description here](https://lh3.googleusercontent.com/YLp-aXY4_5KDswitLovAOqvXhu7r4wOPdjgcIFqVEA_MybzhQqIAUEzA9OIMKAmIYXq-iejOm4TJew)
	
Input = raw counts from htseq-count or featureCounts. Can also use QoRTs Counts - see [page 17](http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf)



```bash
mkdir -p edgeR

GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf

#go to htseq directory & create a mapping file to go from ENSG IDs (which is htseq output) to Symbols:
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $GTF | sort | uniq > ENSG_ID2Name.txt
head ENSG_ID2Name.txt

#determine the number of unique Ensembl Gene IDs & Symbols
## Take the 1st field (column), sort into ASCII order, remove duplicates, print lines + word + byte count
cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc
cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc
## number of occurances of each gene, reverse sort
cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head
```
R script:
```r
Load dataset as above
#################
# Running edgeR #
#################

# load edgeR
library('edgeR')
# make class labels
class <- factor( c( rep("VCP",3), rep("CTRL",3) ))
# Get common gene names
genes=rownames(rawdata)
gene_names=mapping[genes,1]
# Make DGEList object
y <- DGEList(counts=rawdata, genes=genes, group=class)
nrow(y)
# TMM Normalization
y <- calcNormFactors(y)
# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
# Differential expression test
et <- exactTest(y)
# Print top genes
topTags(et)
# Print number of up/down significant genes at FDR = 0.05 significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]
# Output DE genes
# Matrix of significantly DE genes
mat <- cbind(
	genes,gene_names,
	sprintf('%0.3f',log10(et$table$PValue)),
	sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]
colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")
# Order by log fold change
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]
# Save table
write.table(mat, file="edgeR_DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
#To exit R type the following
quit(save="no")
```

Examine the Differentially expressed genes:
```bash
# Look at sigDE genes:
cat edgeR_DE_genes.txt

# pull out the gene IDs
cut -f 2 edgeR_DE_genes.txt | sort  > htseq_edgeR_DE_gene_symbols.txt
```

# Ballgown

https://www.bioconductor.org/packages/release/bioc/html/ballgown.html
Designed to work with StringTie count analysis.

First create a file that lists the expression files, then view that file, then start an R session to examine these results:

```bash
printf "\"ids\",\"type\",\"path\"\n\"UHR_Rep1\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep1\"\n\"UHR_Rep2\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep2\"\n\"UHR_Rep3\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep3\"\n\"HBR_Rep1\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep1\"\n\"HBR_Rep2\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep2\"\n\"HBR_Rep3\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep3\"\n" > UHR_vs_HBR.csv
cat UHR_vs_HBR.csv

R
```
Run the ballgown.R script https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part1_ballgown.R

Output = TSV file. Examine the output:
```bash
# How many passed filter in UHR or HBR?
grep -v feature UHR_vs_HBR_gene_results_filtered.tsv | wc -l

# How many differentially expressed genes were found on this chromosome (p-value < 0.05)?
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | wc -l

# Display the top 20 DE genes. Look at some of those genes in IGV - do they make sense?
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -rnk 3 | head -n 20 #Higher abundance in UHR
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -nk 3 | head -n 20 #Higher abundance in HBR

# Save all genes with P<0.05 to a new file.
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > DE_genes.txt
head DE_genes.txt
```










<!--stackedit_data:
eyJoaXN0b3J5IjpbMTM2MTQxMzgwOSwxMzIzNDM1NjIsMTk2NT
k4NTcyMywxOTMyMzY1Nzg3LC01OTYzMTIxNTUsNjQwNTE3Mzk4
LDYwNzQ0MDgwNywxMTU5Mzc2NzY1LC0xODMxNzA1MjgwLDE3OT
Y1MjY2MjAsMjc2NTM5MjYsNzE4MTIyODUsLTEyNTIwMDIyMzgs
LTE1OTUwNzQxMzYsMjA4NzE1NTEyNywxMjU3MjcyNjIxLC0yMj
UxMTI5MDQsODA0NDM2MDUsNzQ5NjUxNDkzLC0yMTkzNzI0MzZd
fQ==
-->