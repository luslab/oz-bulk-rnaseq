> # Differential Gene Expression Analysis

This chapter covers downstream interpretation of expression & differential estimates. To compare two conditions look at the fraction of transcripts assigned to a specific gene over the total number of reads (total read number differs drastically between samples). The number of sequenced reads depends on:
1. expression level
2. read length
3. sequencing depth
4. expression of all other genes in the sample

Normalisation & Log Transforming Read Counts is performed to ensure that systematic effects not related to the biological differences between samples are removed. Methods of normalising:
- total count
- Counts/million
- **DSeq size factor** using R
- TMM (trimmed mean of M values)
- upper quartile

## Comparison Types

1.  _Within-sample_  comparisons: compare the expression of genes within the same experiment eg in this experiment does gene A express at a higher level than gene B?
2.  _Between-sample_  comparisons: compare the expression of a gene between experimental conditions aka pairwise comparison e.g. has the gene expression for gene A  changed across different experimental conditions? if yes, then differentially expressed (DE).

## Goals of DE analysis
The 2 tasks of DE are to:
1. Estimate the magnitude **fold change** of differential expression in read counts between different conditions, accounting for differences in sequencing depth & variability
2. Estimate the **significance** of the difference, accounting for multiple testing

Null hypothesis = the mean read counts of genes in different samples are equal between different conditions.

Need to estimate the **mean** and **dispersion** from the read counts:
- precision depends on the number & variation of replicates (i.e. the statistical power). 
- For RNA-seq typically there is only 2-3 replicates creating poor precision. There are tools to compensate for this by looking across genes with similar expression to reduce a given genes variance towards the regressed values.

## Tools for DE analysis

DESeq2, edgeR (better for false positives, less conservative, recommended if <12 replicates), ballgown
Ccuffdiff is slow, cant support multifactored experiments, can detect differential isoforms, high false positives.

![enter image description here](https://lh3.googleusercontent.com/LVvCl3GXhNzUx5lyTrHsr0z_ZmI0nb51TBiY1-53VifMuYW8HR9-X54sfLwoH5gFyqahHOm8_QaWhg "Comparison of DGE programs")

DESeq2 and edgeR are similar. Stick to one.

# DESeq2

DESeq takes counts --> read them into R --> normalise for sequencing depth differences

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
ZZ

The results.txt file describes changes between the 2 conditions e.g.
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
#Sort by gene ID select only columns foldchange and log2FoldChange. The results.txt file is already sorted according to padj
cat results.txt | sort | cut -f 1,5,6 > table

#How many genes are significantly differentially expressed (i.e. padj < 0.05 in column 8)?
cat results.txt | awk ' $8 < 0.05 { print $0 }' > diffgenes.txt

#How many differentially expressed genes do we have?
cat diffgenes.txt | wc -l
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



## DSeq2 script

```R
# load magrittr in R & DSeq2
library(magrittr)
library(DESeq2)
# get table of read counts
read.counts = read.table("featureCounts_results.txt", header = TRUE) 
# store gene IDs as row.names
row.names(read.counts) = readcounts$Geneid
# exclude columns that dont have read counts
readcounts = readcounts[ , -c(1:6)]
# assign the sample names
orig_names = names(readcounts)
names(readcounts) = gsub(".*(WT|SNF2)(_[0-9]+).*", "\\1\\2 ", orig_names)
# check data
str(readcounts)
head(readcounts)
# make a data frame for rows (samples)
sample_info = data.frame(condition = gsub("_[0 -9]+", "", names(readcounts)), row.names = names(readcounts))
# Generate the DSeqDataSet - run the DGE analysis
DESeq.ds = DESeqDataSetFromMatrix(countData = readcounts, colData = sample_info, design = ~ condition)
# Check and test dataset
colData(DESeq.ds) %>% head
assay(DESeq.ds) %>% head
rowRanges(DESeq.ds) %>% head
counts(DESeq.ds) %>% str
DESeq.ds = DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]
#remove genes without any counts
colSums(counts(DESeq.ds))` # should be the same as `colSums(readcounts)
# DSeq default for normalising for differences in sequencing depths is `estimateSizeFactors` calculate the size factor and add it to the data set:
DESeq.ds = estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)
# counts ()` allows you to immediately retrieve the normalized read counts
counts.sf_normalized = counts(DESeq.ds, normalized = TRUE)
# Log Transformation of Sequencing Depth Normalised read counts. Most downstream analyses work best on log scales of read counts. Usually *log2* but occasionally *log10*. Log2 transform read counts: 
log.norm.counts = log2(counts.sf_normalized + 1)
#use a pseudocount of 1
```

# edgeR

edgeR is a bioconductor package designed for DE of raw counts
Input = raw counts from htseq-count or featureCounts


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
#######################
# Load Data into R #
#######################

#Set working directory in R where output will go:
working_dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/htseq"
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
eyJoaXN0b3J5IjpbNTAwODgyMzkwLC0xNDkzNzAwNTcxLDE5MT
gxNDA2NTcsLTQ5NzE4NTQxMywyMDIwODg2NzQ4LDkyMDMwNTQ1
NCwyMDM5NzAyODY2LC0xNjQxMTQ1MDEyLDExMjgzODI1MjAsLT
E1MTMzODYzNTUsMTUwNzEzODgwOCwxMjYzOTYxNTY0LDE4Njc3
NTMyMDYsLTIwMzY0MzI3MDcsLTE5MjcwMTAzMjgsLTM1OTYwNT
YzNiwtMTE1NzAwMTA3OCwxNTUyMTcyNTU3LDkxMDE4NDIzMywt
MjEyODIyMzkyNV19
-->