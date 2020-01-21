
> # Visualise Differential Expression

There are many scripts available specific to the tool you used for DE analysis (edgeR, DESeq, ballgown, limma etc). 

There are 2 parts to DE visualisation:
1. Relationships between Samples: Transformations of counts for making plots.
2. Statistical testing for differences attributable to disease & control. 

Check basic global patterns are met:
- do replicates have similar expression patterns 
- do experimental conditions have differences
- For gene which you have prior knowledge about, you should check to see if they behaved as expected - a knockout gene should be very strongly down-regulated in the DGE analysis.
- Map the ORF identifiers from the read count matrix to the gene name --> retrieve the rlog transformed read counts & log2 fold changes.

https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Supplementary_R.R
https://github.com/griffithlab/rnaseq_tutorial/tree/master/scripts

# IGV

Firstly, visualise the most significantly DE genes in IGV
1. Run IGV on local computer and mount CAMP. [Set Java 8 as default](https://stackoverflow.com/questions/46513639/how-to-downgrade-java-from-9-to-8-on-a-macos-eclipse-is-not-running-with-java-9) since IGV doesnt work with Java 10
```bash
# On local terminal move to IGV in bin
cd ~/bin/IGV_2.4.14/lib
# run IGV via command line on local terminal
java -Xmx750m -jar igv.jar
```
2. Set reference genome to Human (hg38) top left box.
3. Click File load from file > click Desktop > mount CAMP locally > click relevant BAM & BAI files or BigWig (.bw) files (can load multiple at once).
4. Search most significantly DE genes reported from DE analysis output in IGV: 
Take ENSG ID > type into Ensembl / Google > Type Gene Name directly in IGV search box



# Visualisation Resources
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#exploratory-analysis-and-visualization

https://github.com/griffithlab/rnaseq_tutorial/wiki/DE-Visualization

http://bioconductor.org/packages/release/bioc/vignettes/ReportingTools/inst/doc/rnaseqAnalysis.pdf

# Results Visualisation Packages
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#rich-visualization-and-reporting-of-results

Start with using an automated visualisation package as it autocreates key plots which are interactive in HTML format. These all work with DESeq2

 - **Glimma.**  Interactive visualization of DESeq2 output, including MA-plots (also called MD-plot) can be generated using the  [Glimma](http://bioconductor.org/packages/Glimma)  package. See the manual page for  _glMDPlot.DESeqResults_.
 - **ReportingTools.**  An HTML report of the results with plots and sortable/filterable columns can be generated using the  [ReportingTools](http://bioconductor.org/packages/ReportingTools)  package on a _DESeqDataSet_ that has been processed by the _DESeq_ function. For a code example, see the _RNA-seq differential expression_ vignette at the [ReportingTools](http://bioconductor.org/packages/ReportingTools) page, or the manual page for the _publish_ method for the _DESeqDataSet_ class.
 - **regionReport.**  An HTML and PDF summary of the results with plots can also be generated using the  [regionReport](http://bioconductor.org/packages/regionReport)  package. The  _DESeq2Report_  function should be run on a  _DESeqDataSet_  that has been processed by the  _DESeq_  function. For more details see the manual page for  _DESeq2Report_  and an example vignette in the  [regionReport](http://bioconductor.org/packages/regionReport)  package.
 - **pcaExplorer.**  Interactive visualization of DESeq2 output, including PCA plots, boxplots of counts and other useful summaries can be generated using the  [pcaExplorer](http://bioconductor.org/packages/pcaExplorer)  package. See the  *Launching the application* section of the package vignette.

```r
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
```

## regionReport
http://leekgroup.github.io/regionReportSupp/DESeq2.html
[http://bioconductor.org/packages/release/bioc/vignettes/regionReport/inst/doc/regionReport.html](http://bioconductor.org/packages/release/bioc/vignettes/regionReport/inst/doc/regionReport.html)
```r
library('ggplot2')
library('regionReport')
report <- DESeq2Report(dds, project = 'DESeq2 HTML report', 
	intgroup = c('condition'), outdir = 'regionReport', 
	output = 'index', theme = theme_bw())
```

Now create plots manually. SVD analysis.Rmd script runs the visualisaiton manually.

# Transformations
Statistical exploration of multidimension data eg clustering & PCA, work best when data has the same range of variance across the range of means (homoskedastic). With RNA-seq raw counts the variance grows with the mean. DESeq2 uses **rlog** (regularise logarithm) & **vst** (variance stabilising transformation) to transform count data stabilising variance across the mean. 

Continue in R from DE analysis using QoRTS > DESeq2:
```R
#Set working directory where counts expression files exist
working_dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2"
setwd(working_dir)

res <- read.qc.results.data("/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/", 
                            decoder.files = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.byUID.txt", 
                            calc.DESeq2 = TRUE, calc.edgeR = TRUE); 

### Run DESeq2 DE analysis
decoder.bySample <- read.table("/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.bySample.txt", 
                               header=T,stringsAsFactors=F); 
directory <- "/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC"; 
sampleFiles <- paste0(decoder.bySample$qc.data.dir, "/QC.geneCounts.formatted.for.DESeq.txt.gz" ); 
sampleCondition <- decoder.bySample$group.ID; 
sampleName <- decoder.bySample$sample.ID; 
sampleTable <- data.frame(sampleName = sampleName, 
                          fileName = sampleFiles, 
                          condition = sampleCondition); 
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                  directory = directory, 
                                  design = ~ condition); 
dds <- DESeq(dds);
res <- results(dds);
resLFC <- lfcShrink(dds, coef="condition_VCP_vs_CTRL", type="apeglm")

### Data Transformation
# extract transformed values using vst (rapid) or rlog (slower)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
# extract the matrix of normalised values
head(assay(vsd), 6)
head(assay(rld), 6)
# print column metadata in DESeqDataSet
colData(vsd)
colData(rld)
# Visually assess effect of transformation on variance. Plot SD of transformed data vs mean
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

```
# PCA Plot
[https://www.youtube.com/watch?v=HMOI_lkzW08](https://www.youtube.com/watch?v=HMOI_lkzW08)
[https://www.youtube.com/watch?v=FgakZw6K1QQ](https://www.youtube.com/watch?v=FgakZw6K1QQ)

Visualise sample-sample distances with the Principal Components Analysis (PCA). Data points represents samples in 2D spread in 2 directions: x-axis = PC1; y-axis = PC2. % of total variance is printed in axis label. These dont add to 100% as there are more dimension containing remaining variance. 

Assess if samples have greater variance between experimental and control conditions than between replicates. Result is principal components representing directions along which the variation in the originial multi-dimensional data is maximal, so that a few components (dimensions) can be used to represent thousands of mRNA data points. Can visually represent variation of gene expression for different samples in a simple xy plot (instead of plotting thousands of genes per sample). Usually only the top 2 principal components (explaining the majority of the data variability) are displayed.

![enter image description here](https://onlinecourses.science.psu.edu/stat857/sites/onlinecourses.science.psu.edu.stat857/files/lesson05/PCA_plot/index.gif)

```r
#quick PCA
plotPCA(vsd, "condition")
plotPCA(rld, "condition")
plotPCA(ntd, "condition")

# using ggplot2 transformed with VST (can change to rlog or ntd)
data <- plotPCA(vsd, intgroup = c( "condition", "sizeFactor"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=name)) + geom_point(size=3) + 
	ggtitle("PCA Plot VST transformed counts") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
```

# MDS Plot

Multidimensional scaling (MDS) plot. Similar to PCA as reduces dimensions. Calculate pairwise distances between samples then crease a 2D representation of these distances.
```r
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)), decoder.bySample)
head(mds)
ggplot(mds, aes(X1,X2,color=condition,shape=sample.ID)) + geom_point(size=3)+ 
  ggtitle("MDS Plot VST transformed counts")
```

# Single gene Plot Counts
It can be useful to examine the counts of reads for a single gene across the groups. The function for making this plot is  _plotCounts_, which normalizes counts by sequencing depth and adds a pseudocount of 1/2 to allow for log scale plotting. 

```r
#quickly visualise
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, topGene, "condition")

#with ggplot
d <- plotCounts(dds, gene=rownames(res)[which.min(res$padj)], intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```

# MA plot
An MA plot provides a global overview of an experiment with a 2 group comparison. It shows the relationship between expression change in 2 different conditions.  y-axis = log2 fold change for comparison. x-axis = average counts normalised by size factor. Each gene is a dot. Genes with p-adj <0.05 are shown in red.

![enter image description here](https://lh3.googleusercontent.com/PdRsM9aHl3MTvEMKCYjYKQysVZ9MKxk943_XZ_JLLtAH0jTgZXKP2XotWhetjvghPqGDdwn0ULGRBw "Histogram & MA plot")
```r
DESeq2::plotMA(res, alpha=0.05, main = "MA plot of VCP vs. CTRL expression", ylim=c(-6,6))

# Plot for shrunken log2 fold changes (removes noise from log2 fold changes from low count genes)
DESeq2::plotMA(resLFC, alpha=0.05, main = "MA plot of VCP vs. CTRL expression VCP vs. CTRL", ylim=c(-6,6))
```

# Heatmap
https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#heatmap-of-the-most-significant-genes

Data quality testing is essential early in the analysis. Remove poor data that suffers from anormality. Show expression values across individual samples. many R functions to do this:
`aheartmap( )`
`gplots::heatmap.2( )`
`pheatmap::pheatmap( )`

Construct a heatmap to explore the count matrix:
```r
library(pheatmap)
### Heatmap of most significantly DE genes
mat <- assay(vsd)[ head(order(res$padj),30), ]
mat <- mat - rowMeans(mat)
df <- cbind(as.data.frame(colData(vsd), decoder.bySample$sample.ID))
pheatmap(mat, annotation_col=df, main = "Heatmap: Most significantly DE genes")
```

Biostars code to generate a clustered heatmap: 
`curl -O http://data.biostarhandbook.com/rnaseq/code/draw-heatmap.r`
This generates a PDF output using:
`cat normalise-matrix-deseq.txt | Rscript draw-heatmap.r > clustered-heatmap.pdf`
Each column referrs to a sample. Red refers to upregulated genes & green downregulated.

![enter image description here](http://bioinfo.cipf.es/babelomicstutorial/_media/images:differential_expression_example:heatmap.png)


## Hierarchical clustering
- separate samples in an **unsupervised** fashion to see if samples of different conditions differ more than replicates within the same condition.
- pairwise compairsons of individual samples, grouped into "neighbourhoods" of similar samples.
- hierachical clustering analyses require decisions on:
	- how the disimilarity between pairs should be calculated?
	- how the disimilarity be used for clustering?
		- Pearson correlation coefficient
		- Euclidean distance (distance between two vectors) - calculate the distance using `linkage` function (complete, average, or single intercluster distance - complete intercluster distance is best and single IC distance is worst).

Can visualise DE genes as a Dendogram or Heatmap

## Dendogram
Genes are sorted by adjusted p-value. Colours represent read counts.
- clustures obtained by cutting dendoram at a level where the jump between two nodes is large
- connected components form individual clusters
- clustering algorithms differ and there is no concensus on which is optimal
in R use `cor( )` `as.dist( )` and `hclust( )` to generate a dendogram



## Annotate & export results:
https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#annotating-and-exporting-results
```r
library("AnnotationDbi")
library("Homo.sapiens")
### Annotation
columns(Homo.sapiens)
# remove decimal string & everything that follows from row.names (ENSEMBL) and call it "tmp": https://www.biostars.org/p/178726/
tmp=gsub("\\..*","",row.names(res))
# use mapIds to add columns to results table. Use Ensebml nomenclature. 
res$symbol <- mapIds(Homo.sapiens,
                     keys=tmp,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
# ignore the error: 'select()' returned 1:many mapping between keys and columns
res$entrez <- mapIds(Homo.sapiens,
                     keys=tmp,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$symbol <- mapIds(Homo.sapiens,
                     keys=tmp,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$padj),]
head(resOrdered)

#view the dataframe with ENSEMBL symbol external gene IDs
resOrderedDF <- as.data.frame(resOrdered)
col_symbol <- grep("symbol", names(resOrderedDF))
resOrderedDF <- resOrderedDF[, c(col_symbol, (1:ncol(resOrderedDF))[-col_symbol])]
names(resOrderedDF)
head(resOrderedDF)
View(resOrderedDF)
# can export top 100 DE genes to CSV file
resOrderedDF_top100 <- as.data.frame(resOrdered)[seq_len(100),]
write.csv(resOrderedDF_top100, file="DESeq2_results.csv")
```

# Number of transcripts per gene

Use output from DE analysis

```R
#### Plot #1 - the number of transcripts per gene.
#Many genes will have only 1 transcript, some genes will have several transcripts
#Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol)
#Then use the 'hist' command to create a histogram of these counts
#How many genes have 1 transcript? More than one transcript? What is the maximum number of transcripts for a single gene?
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
```

# Transcript sizes as Histogram

```r
#### Plot #2 - the distribution of transcript sizes as a histogram
#In this analysis we supplied StringTie with transcript models so the lengths will be those of known transcripts
#However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled
#If we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts
full_table <- texpr(bg , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

#### Summarize FPKM values for all 6 replicates
#What are the minimum and maximum FPKM values for a particular library?
min(gene_expression[,"FPKM.UHR_Rep1"])
max(gene_expression[,"FPKM.UHR_Rep1"])
#Set the minimum non-zero FPKM values for use later.
#Do this by grabbing a copy of all data values, coverting 0's to NA, and calculating the minimum or all non NA values
#zz = fpkm_matrix[,data_columns]
#zz[zz==0] = NA
#min_nonzero = min(zz, na.rm=TRUE)
#min_nonzero
#Alternatively just set min value to 1
min_nonzero=1
# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:6)
short_names=c("UHR_1","UHR_2","UHR_3","HBR_1","HBR_2","HBR_3")
```

# DE Histogram
```r

#### Plot #9 - View the distribution of differential expression values as a histogram
#Display only those that are significant according to Ballgown
sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) UHR vs HBR", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

#### Plot #10 - Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed
gene_expression[,"UHR"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"HBR"]=apply(gene_expression[,c(4:6)], 1, mean)
x=log2(gene_expression[,"UHR"]+min_nonzero)
y=log2(gene_expression[,"HBR"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="UHR FPKM (log2)", ylab="HBR FPKM (log2)", main="UHR vs HBR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"qval"])[1:25]
text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)

#### Write a simple table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

#Order the output by or p-value and then break ties using fold-change
o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]

#You can open the file "SigDE.txt" in Excel, Calc, etc.
#It should have been written to the current working directory that you set at the beginning of the R tutorial
dir()
```

# Range of FPKM
```r
#### Plot #3 - View the range of values and general distribution of FPKM values for all 4 libraries
#Create boxplots for this purpose
#Display on a log2 scale and add the minimum non-zero value to avoid log2(0)
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 6 libraries")
#Note that the bold horizontal line on each boxplot is the median
```
# Replicates
```r
#### Plot #4 - plot a pair of replicates to assess reproducibility of technical replicates
#Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0)
x = gene_expression[,"FPKM.UHR_Rep1"]
y = gene_expression[,"FPKM.UHR_Rep2"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (UHR, Replicate 1)", ylab="FPKM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates")

#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)

#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
```
# Density scatter plot
```r

#### Plot #5 - Scatter plots with a large number of data points can be misleading ... regenerate this figure as a density scatter plot
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (UHR, Replicate 1)", ylab="FPKM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)

#### Plot all sets of replicates on a single plot
#Create an function that generates an R plot. This function will take as input the two libraries to be compared and a plot name and color
plotCor = function(lib1, lib2, name, color){
x=gene_expression[,lib1]
y=gene_expression[,lib2]
zero_count = length(which(x==0)) + length(which(y==0))
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col=color, cex=0.25, xlab=lib1, ylab=lib2, main=name)
abline(a=0,b=1)
rs=cor(x,y, method="pearson")^2
legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

#Open a plotting page with room for two plots on one page
par(mfrow=c(1,2))

#Plot #6 - Now make a call to our custom function created above, once for each library comparison
plotCor("FPKM.UHR_Rep1", "FPKM.HBR_Rep1", "UHR_1 vs HBR_1", "tomato2")
plotCor("FPKM.UHR_Rep2", "FPKM.HBR_Rep2", "UHR_2 vs HBR_2", "royalblue2")

##### One problem with these plots is that there are so many data points on top of each other, that information is being lost

#Regenerate these plots using a density scatter plot
plotCor2 = function(lib1, lib2, name, color){
x=gene_expression[,lib1]
y=gene_expression[,lib2]
zero_count = length(which(x==0)) + length(which(y==0))
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
abline(a=0,b=1)
rs=cor(x,y, method="pearson")^2
legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

#### Plot #7 - Now make a call to our custom function created above, once for each library comparison
par(mfrow=c(1,2))
plotCor2("FPKM.UHR_Rep1", "FPKM.HBR_Rep1", "UHR_1 vs HBR_1", "tomato2")
plotCor2("FPKM.UHR_Rep2", "FPKM.HBR_Rep2", "UHR_2 vs HBR_2", "royalblue2")

#### Compare the correlation 'distance' between all replicates
#Do we see the expected pattern for all eight libraries (i.e. replicates most similar, then tumor vs. normal)?

#Calculate the FPKM sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")

#Print out these correlation values
r
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTgxNDY1ODUxMCwtMTk0NTIwNDYxNl19
-->