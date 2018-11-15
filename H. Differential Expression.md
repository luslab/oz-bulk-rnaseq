


> # Differential Gene Expression Analysis

# Comparisons between conditions

To compare two conditions look at the fraction of transcripts assigned to a specific gene over the total number of reads (total read number differs drastically between samples). The number of sequenced reads depends on:
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
The 2 tasks of DGE are to:
1. Estimate the magnitude **fold change** of differential expression in read counts between different conditions, accounting for differences in sequencing depth & variability
2. Estimate the **significance** of the difference, accounting for multiple testing

Null hypothesis = the mean read counts of genes in different samples are equal between different conditions.

Need to estimate the **mean** and **dispersion** from the read counts:
- precision depends on the number & variation of replicates (i.e. the statistical power). 
- For RNA-seq typically there is only 2-3 replicates creating poor precision. There are tools to compensate for this by looking across genes with similar expression to reduce a given genes variance towards the regressed values.

## Tools for DE analysis

There are 2 approaches:
1. Rapid approach: [Kallisto](https://pachterlab.github.io/kallisto/)  which is a really user-friendly algo which extract both gene and transcript level gene expression directly from fastq files using the raw fastq files. Then use  [Sleuth](https://pachterlab.github.io/sleuth/)  (also developed by Pachter lab) to perform differential gene and transcript expression analysis.
2. Detailed approach: DESeq, DESeq2, edgeR

`edgeR` (better for false positives, less conservative, recommended if <12 replicates)
`DESeq`
`DESeq2` sample wise size factor
`limma-voom`
`cuffdiff`slow, cant support multifactored experiments, can detect differential isoforms, high false positives

![enter image description here](https://lh3.googleusercontent.com/LVvCl3GXhNzUx5lyTrHsr0z_ZmI0nb51TBiY1-53VifMuYW8HR9-X54sfLwoH5gFyqahHOm8_QaWhg "Comparison of DGE programs")

# Rapid Approach: Kallisto - Sleuth pipeline
Author = [Lior Patcher](https://en.wikipedia.org/wiki/Lior_Pachter)

Kallisto quantifies transcript abundances. It pseudoaligns reads against a transcriptome (not genome). Simple count-based approaches underperform when determining transcript level counts as they disregard reads that overlap with more than one gene. If the genomic feature becomes a transcript rather than a gene it keeps many reads that would have been discarded.

Method of pseudoalignment: for each read the program aims to identify the target that it originates from using k-mers. By ignoring exactly where in the genome a read originates from it is much faster than normal alignment. This approach does not generate a BAM file (alignment file) but instead produce a measure of how many reads indicate the presence of each transcript. These use a **deBruikin graph** to assign reads to an isoform if they are compatible with that transcript structure.
![enter image description here](https://www.frontiersin.org/files/Articles/169488/fgene-06-00361-r2/image_m/fgene-06-00361-g002.jpg)

Schema of a simple deBruijn graph-based transcript assembly:
- Read sequences are split into all subsequence k-mers (here: of length 5) from the reads.
- A deBruijn graph is constructed using unique k-mers as the nodes and overlapping k-mers connected by edges (a k-mer shifted by one base overlaps another k-mer by k􀀀1 bases).
- The transcripts are assembled by traversing the two paths in the graph

Although much faster than alignment-counting routines it **cant detect novel isoforms**. However, instead of direct isoform quantification, you can glean more accurate answers from alternative approaches, e.g., quantification of exons (Anders et al., 2012) or estimates of alternative splicing events such as exon skipping, intron retention etc. (e.g., MISO [Katz et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037023/), rMATS [Shen et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280593/)).

The main limitations to assigning reads to transcripts are:
- annotation transcripts are inconsistent
- many isoforms with very different lengths
- anti-sense and overlapping transcripts of different genes

**Tools**:
Sailfish and more updated version Salmon
[Kallisto](https://www.nature.com/articles/nbt.3519)

## Kallisto Workflow
ml kallisto

```bash
#set changable elements
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.fa
IDX=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.cdna.fa.idx
SAMPLE=CTRL_1
R1=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_depleted/${SAMPLE}.fq
##if you have paired-end data then set the 2nd read files for input
R2=PATH_TO_FASTQ_reverse_strand.fq

#build kallisto index
sbatch -N 1 -c 8 --mem 40 --wrap="kallisto index -i $IDX $REF"

#create directory for kallisto data 
mkdir -p kallisto
#set subdirectory for output called out
OUTDIR=out
#run kallisto quantification with quant command. -o sets output directory. -b specifies the bootstap sample number.
##paired-end mode
kallisto quant -i $IDX -o $OUTDIR -b 100 $R1 $R2
##single-end mode (set fragment mean length & standard deviation - Illumina generates fragments 180-200bp - acurately determine this from a library quantification with an instrument such as an Agilent Bioanalyzer)
kallisto quant -i $IDX -o $SAMPLE -b 100 --single -l 187 -s 70 $R1

#move the output to kallisto results folder
mv $SAMPLE /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto
```
You can run this for multiple samples as a For Loop (example [here](https://www.biostarhandbook.com/rnaseq/rnaseq-griffith-kallisto.html))
```bash
# Create output folder
mkdir -p kallisto

# Exit this script on any error.
set -euo pipefail

# Set Reference transcriptome (not genome).
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.fa
# Set index to build
IDX=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.cdna.fa.idx

# Build kallisto index
sbatch -N 1 -c 8 --mem=40GB --wrap="kallisto index -i $IDX  $REF"

for SAMPLE in VCP CTRL;
do
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files.
        R=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_depleted/${SAMPLE}_${REPLICATE}.fq
        # The kallisto output directory.
        OUTDIR=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/${SAMPLE}_${REPLICATE}
        # Run kallisto quantification in single-end mode.
        kallisto quant --single -l 187 -s 70 -i $IDX -o $OUTDIR -b 100 $R

        # Copy the abundance file to a proper name - i.e. remove the long path name so that it only contains information on the sample e.g. VCP_3.tsv
        cp $OUTDIR/abundance.tsv $OUTDIR.counts.tsv
    done
done

#concatenate all counts into 1 file that can be used for DESeq DE analysis.
paste /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/VCP_*.tsv  /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/CTRL_*.tsv | cut -f 1,4,9,14,19,24,29 > counts.txt
```
### Kallisto Output

The main output of Kallisto is the abundance.txt file with columns:
```
target_id   length eff_length est_counts 	tpm
ERCC-00002  1061   891.059      18946      243099
ERCC-00003  1023   853.059      1452       19460.7
ERCC-00004  523    353.059      455        14734.5
ERCC-00009  984	   814.059      319        4480.29
```
eff_length = scales the transcript length by fragment length distribution . (transcript length - mean fragment length + 1)
est_counts = transcript abundances
tpm = Transcripts Per Million

# Detailed Approach: DESeq2 > Visualisation

DESeq( ) function wraps around the following 3 functions:
`DESeq.ds = estimateSizeFactors(DESeq.ds)` #sequencing depth normalisation
`DESeq.ds = estimateDispersions(DESeq.ds)` #gene wise dispersion estimates
`DESeq.ds = nbinomWaldTest(DESeq.ds)` #fits a negative binomial GLM & applies Wald stats to each gene

`DGE. results = results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)` #results function allows extraction of base means across sample, SEs, and basic stats for each gene.
`summary(DGE.results)`

The DESeqResult behaves as a data frame:
head (DGE. results )
table (DGE. results $ padj < 0.05)
rownames ( subset (DGE. results , padj < 0.05) )

 Takes the `featureCounts` (raw read counts) --> read them into R --> normalise for sequencing depth differences
 
 `DSeqDataSet` contains all the information in R. 
`rowRanges( )` rows = variables (genes, transcripts, exons)  - has info about genes (chr, start, end, strand, ID)
`rownames` is the unique sample names
`colData( )` columns = samples 
`assay( )` stores count data with genes and samples. similar to `countData`

## DSeq2 process

1. load magrittr in R & DSeq2:
`library(magrittr)`
`library(DESeq2)`
2. get table of read counts:
`read.counts = read.table("featureCounts_results.txt", header = TRUE)` 
3. store gene IDs as row.names:
`row.names(read.counts) = readcounts$Geneid`
4. exclude columns that dont have read counts:
`readcounts = readcounts[ , -c(1:6)]`
5. assign the sample names:
`orig_names = names(readcounts)`
`names(readcounts) = gsub(".*(WT|SNF2)(_[0-9]+).*", "\\1\\2 ", orig_names)`
6. check data:
`str(readcounts)`
`head(readcounts)`
7. make a data frame for rows (samples)
`sample_info = data.frame(condition = gsub("_[0 -9]+", "", names(readcounts)), row.names = names(readcounts))`
8. Generate the DSeqDataSet - run the DGE analysis
`DESeq.ds = DESeqDataSetFromMatrix(countData = readcounts, colData = sample_info, design = ~ condition)`
9. Check and test dataset
`colData(DESeq.ds) %>% head`
`assay(DESeq.ds) %>% head`
`rowRanges(DESeq.ds) %>% head` 
`counts(DESeq.ds) %>% str`
`DESeq.ds = DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]` #remove genes without any counts
`colSums(counts(DESeq.ds))` # should be the same as `colSums(readcounts)`

DSeq default for normalising for differences in sequencing depths is `estimateSizeFactors`
calculate the size factor and add it to the data set:
`DESeq.ds = estimateSizeFactors(DESeq.ds)`
`sizeFactors(DESeq.ds)`
`counts ()` allows you to immediately retrieve the normalized read counts:
`counts.sf_normalized = counts(DESeq.ds, normalized = TRUE)`

10. Log Transformation of Sequencing Depth Normalised read counts

Most downstream analyses work best on log scales of read counts. Usually *log2* but occasionally *log10*.
Log2 transform read counts: `log.norm.counts = log2(counts.sf_normalized + 1)` #use a pseudocount of 1

## Biostars DESeq Script
```bash
ml R

#download the DESeq biostars script
curl -O http://data.biostarhandbook.com/rnaseq/code/deseq1.r
curl -O http://data.biostarhandbook.com/rnaseq/code/deseq2.r

# Using featureCounts output (counts.txt) print out only columns representing Gene ID & sample abundances (ie remove intermediate columns - Chr, Start, End, Strand & length)
cat counts.txt | cut -f 1,7-14 > sample_counts.txt

#pass simple_counts.txt through the script specifying the design of the experiment 
## in this case = 3 x 3 (3 cases, 3 controls)
cat sample_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.txt
```
The results_deseq1.txt file describes changes between the 2 conditions e.g.
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

# Visualise Differential Expression

## IGV

Firstly, visualise the most significantly DE genes in IGV
1. On local terminal `cd ~/bin/IGV_2.4.14/lib` & run IGV via command line on local terminal: `java -Xmx750m -jar igv.jar`
2. Set reference genome to Human (hg38) top left box.
3. Click File load from file > click Desktop > mount CAMP locally > click relevant BAM & BAI files (can load multiple at once).


## Explore Read Count Patterns

Check that basic global patterns are met:
- do replicates have similar expression patterns; 
- do experimental conditions have differences
- For gene which you have prior knowledge about, you should check to see if they behaved as expected - a knockout gene should be very strongly down-regulated in the DGE analysis.
- Map the ORF identifiers from the read count matrix to the gene name --> retreive the rlog transformed read counts & log2 fold changes.

Assess RNA-seq expression patterns with:
1. Pairwise Correlation
2. Hierarchical Clustering
3. Principal Component Analysis

## Pairwise correlation

- Pearson correlation coefficient, r is used to assess similarity between RNA-seq samples in a pair wise fashion.
- ENCODE guideline advises **>0.9** correlation should achieved for mRNA transcripts.
- in R use `cor( )` function

## Hierarchical clustering
- separate samples in an **unsupervised** fashion to see if samples of different conditions differ more than replicates within the same condition.
- pairwise compairsons of individual samples, grouped into "neighbourhoods" of similar samples.
- hierachical clustering analyses require decisions on:
	- how the disimilarity between pairs should be calculated?
	- how the disimilarity be used for clustering?
		- Pearson correlation coefficient
		- Euclidean distance (distance between two vectors) - calculate the distance using `linkage` function (complete, average, or single intercluster distance - complete intercluster distance is best and single IC distance is worst).

Can visualise DE genes as a Dendogram or Heatmap

### Clustered Heatmap
Very popular method.

- show expression values across individual samples
many R functions to do this:
`aheartmap( )`
`gplots::heatmap.2( )`
`pheatmap::pheatmap( )`

Biostars code to generate a clustered heatmap: 
`curl -O http://data.biostarhandbook.com/rnaseq/code/draw-heatmap.r`
This generates a PDF output using:
`cat normalise-matrix-deseq.txt | Rscript draw-heatmap.r > clustered-heatmap.pdf`
Each column referrs to a sample. Red refers to upregulated genes & green downregulated.

![enter image description here](http://bioinfo.cipf.es/babelomicstutorial/_media/images:differential_expression_example:heatmap.png)

### Dendogram
Genes are sorted by adjusted p-value. Colours represent read counts.
- clustures obtained by cutting dendoram at a level where the jump between two nodes is large
- connected components form individual clusters
- clustering algorithms differ and there is no concensus on which is optimal
in R use `cor( )` `as.dist( )` and `hclust( )` to generate a dendogram




![enter image description here](http://www.sthda.com/english/sthda-upload/figures/cluster-analysis/009c-divisive-hierarchical-clustering-compute-diana-1.png)

### Principal Components Analysis (PCA)

 - Complementary approach to assess if samples have greater variance between experimental and control conditions than between replicates.
 - Aim is to **identify groups of features** (eg genes) that have something in common, such as expression patterns across different samples.
 - Result is principal components representing directions along which the variation in the originial multi-dimensional data is maximal, so that a few components (dimensions) can be used to represent thousands of mRNA data points.
 - Can visually represent variation of gene expression for different samples in a simple xy plot (instead of plotting thousands of genes per sample)/ Usually only the top 2 principal components (explaining the majority of the data variability) are displayed.
	 - identify unexpected patterns - batch effects; outliers
	 - does not identify unknown groupings

in R use `prcomp` function"
`library(DESeq2)`
`library(ggplot2)`
`pc = prcomp(t(rlog.norm.counts))`
`plot(pc$x[ ,1], pc$x[ ,2], col = colData(DESeq.ds)[ ,1], main = "PCA of seq.depth normlised\n and rlog-transformed read counts"`
`P = plotPCA(DESeq.rlog)` # PCA plot using DESeq2 based on ggplot2
`P = P + theme_bw() + ggtitle("Rlog transformed counts")` #plot cosmetics
`print(P)`

![enter image description here](https://onlinecourses.science.psu.edu/stat857/sites/onlinecourses.science.psu.edu.stat857/files/lesson05/PCA_plot/index.gif)

## Explore normalised read counts
`par(mfrow = c(2, 1))` #plot the 2 image on top of each other
`boxplot(counts.sf_normalized, notch = TRUE, main = "untransformed read counts", ylab = "read counts")` #plots the non-logged boxplots
`boxplot(log.norm.counts, notch = TRUE, main = "log2 - transformed read counts", ylab = "log2(read counts)")` #plots the logged boxplots

Plot replicate results in a pairwise manner
`plot(log.norm.counts[ ,1:2], cex =.1, main = "Normalized log2 (read counts)")`

Check data to ensure variable have similar variance (homoskedastic behaviour):
`library(vsn)`
`library(ggplot2)`
`msd_plot = meanSdPlot(log.norm.counts, ranks =FALSE , plot = FALSE)` #plot mean against SD
`msd_plot$gg + ggtitle ("sequencing depth normalized log2 (read counts)") + ylab("standard deviation")`
y-axis shows variance of read counts. Any rise in the best fit line indicates an increase in variance at that read count length (x axis) - if to left = shorter count lengths; to right = long count lengths.

## Histogram

hist (DGE. results $ pvalue, col = " grey ", border = " white ", xlab = "", ylab = "", main = " frequencies of p- values ") # histogram of p-valus

## MA plot
- provides global view of relationship between expression change in different conditions, average expression strength of genes and ability og algorithrm to detect differential expression
	- red dots are significant adjusted p<0.05
`plotMA(DGE.results, alpha = 0.05, main = "WT vs. SNF2 mutants", ylim = c(-4, 4))`

![enter image description here](https://lh3.googleusercontent.com/PdRsM9aHl3MTvEMKCYjYKQysVZ9MKxk943_XZ_JLLtAH0jTgZXKP2XotWhetjvghPqGDdwn0ULGRBw "Histogram & MA plot")

## Variance Shrinkage

`DSeq2` and `edgeR` both offer means to reduce the variance using the dispersion mean trend using the entire dataset as a reference.

Low read counts that are highly variable will be assigned more homogenous read estimates --> variance resembles the majority of the genes and hopefully has a more stable variance

Regularise log-transformed values:
`DSeq.rlog = rlog(DESeq.ds, blind = TRUE)` #can set rlog to FALSE if there are large differences in a large proportion of the genes to avoid overestimating the dispersion
`rlog.norm.counts = assay(DESeq.rlog)`

`msd_plot = meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)` #show data on original scale and dont print plot
`msd_plot$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation)`


**SVD (singular value decomposition) analysis**

-   For doing this you can use the gene-level count table obtained from Kallisto. I wrote everything in R and I can send you some litterature which explains a bit the underlying math and idea. Also happy to speak about it over skype.


 










<!--stackedit_data:
eyJoaXN0b3J5IjpbOTYyNzk1MzM0LDEwNjA5OTgwNjYsLTE0MD
IzNTEzNzQsNzMwMzE4NTkxLDY2MzY5NjkwMywtMzIyMzg2MzU2
LDYxMTU5MjkzMiwxMjM4MjYwODg4LC0xMjUxNDA1MzU1LC0xNT
E5MTExMDk4XX0=
-->