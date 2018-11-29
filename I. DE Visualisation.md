


> # Visualise Differential Expression

There are many scripts available specific to the tool you used for DE analysis (edgeR, DESeq, ballgown etc). 

There are 2 parts to DE visualisation:
1. Relationships between Samples: Transformations of counts for making plots.
2. Statistical testing for differences attributable to disease & control. 

Check that basic global patterns are met:
- do replicates have similar expression patterns; 
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

# Results Visualisation Packages

These all work with DESeq2:

**ReportingTools.**  An HTML report of the results with plots and sortable/filterable columns can be generated using the  [ReportingTools](http://bioconductor.org/packages/ReportingTools)  package on a _DESeqDataSet_ that has been processed by the _DESeq_ function. For a code example, see the _RNA-seq differential expression_ vignette at the [ReportingTools](http://bioconductor.org/packages/ReportingTools) page, or the manual page for the _publish_ method for the _DESeqDataSet_ class.

```r
desReport <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', 
		title = 'RNA-seq analysis of differential expression using DESeq2',
		reportDirectory = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2/reports") 
publish(res,desReport,name="df",countTable=mockRnaSeqData, pvalueCutoff=0.05,
		conditions=conditions,annotation.db="org.Mm.eg.db", 
		expName="deseq",reportDir="./reports", .modifyDF=makeDESeqDF)
finish(desReport)
```

**regionReport.**  An HTML and PDF summary of the results with plots can also be generated using the  [regionReport](http://bioconductor.org/packages/regionReport)  package. The  _DESeq2Report_  function should be run on a  _DESeqDataSet_  that has been processed by the  _DESeq_  function. For more details see the manual page for  _DESeq2Report_  and an example vignette in the  [regionReport](http://bioconductor.org/packages/regionReport)  package.

**Glimma.**  Interactive visualization of DESeq2 output, including MA-plots (also called MD-plot) can be generated using the  [Glimma](http://bioconductor.org/packages/Glimma)  package. See the manual page for  _glMDPlot.DESeqResults_.

**pcaExplorer.**  Interactive visualization of DESeq2 output, including PCA plots, boxplots of counts and other useful summaries can be generated using the  [pcaExplorer](http://bioconductor.org/packages/pcaExplorer)  package. See the  *Launching the application* section of the package vignette.

# Transformations
Statistical exploration of multidimension data eg clustering & PCA, work best when data has the same range of variance across the range of means (homoskedastic). With RNA-seq raw counts the variance grows with the mean. DESeq2 uses **rlog** (regularise logarithm) & **vst** (variance stabilising transformation) to transform count data stabilising variance across the mean. 



# DESeq2 Visualisation
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#exploratory-analysis-and-visualization
https://github.com/griffithlab/rnaseq_tutorial/wiki/DE-Visualization
http://bioconductor.org/packages/release/bioc/vignettes/ReportingTools/inst/doc/rnaseqAnalysis.pdf

Continue in R from DE analysis using QoRTS > DESeq2:
```R
#Load libraries
library(QoRTs)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library("vsn")

#### Import the gene expression data
#Set working directory where results files exist
working_dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/htseq"
setwd(working_dir)
# List the current contents of this directory
dir()
```
## Clustering

Data quality testing is essential early in the analysis. Remove poor data that suffers from anormality.

Construct a heatmap to explore the count matrix
```r
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```





FROM BALLGOWN GRIFFITH ANALYSIS: REMOVE AFTER USING
```r

#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
load('bg.rda')

# View a summary of the ballgown object
bg

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))

#### Working with 'dataframes'
#View the first five rows of data (all columns) in one of the dataframes created
head(gene_expression)

#View the column names
colnames(gene_expression)
#View the row names
row.names(gene_expression)

#Determine the dimensions of the dataframe. 'dim()' will return the number of rows and columns
dim(gene_expression)

#Get the first 3 rows of data and a selection of columns
gene_expression[1:3,c(1:3,6)]

#Do the same thing, but using the column names instead of numbers
gene_expression[1:3, c("FPKM.UHR_Rep1","FPKM.UHR_Rep2","FPKM.UHR_Rep3","FPKM.HBR_Rep3")]

#Assign colors to each. You can specify color by RGB, Hex code, or name

#To get a list of color names:
colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

#View expression values for the transcripts of a particular gene symbol of chromosome 22. e.g. 'TST'
#First determine the rows in the data.frame that match 'TST', aka. ENSG00000128311, then display only those rows of the data.frame
i = row.names(gene_expression) == "ENSG00000128311"
gene_expression[i,]

#What if we want to view values for a list of genes of interest all at once?
#genes_of_interest = c("TST", "MMP11", "LGALS2", "ISX")
genes_of_interest = c("ENSG00000128311","ENSG00000099953","ENSG00000100079","ENSG00000175329")
i = which(row.names(gene_expression) %in% genes_of_interest)
gene_expression[i,]

# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
length(row.names(transcript_gene_table)) #Transcript count
length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count
```

# MA plot
- provides global view of relationship between expression change in different conditions, average expression strength of genes and ability of algorithrm to detect differential expression

```r
# Plot log2 fold changes over the mean of normalised counts for all samples. red dots are significant adjusted p<0.05
plotMA(res, alpha=0.05, main = "VCP vs. CTRL", ylim=c(-4,4))

# Plot for shrunken log2 fold changes (removes noise from log2 fold changes from low count genes)
plotMA(resLFC, alpha=0.05, main = "VCP vs. CTRL", ylim=c(-4,4))


```
![enter image description here](https://lh3.googleusercontent.com/PdRsM9aHl3MTvEMKCYjYKQysVZ9MKxk943_XZ_JLLtAH0jTgZXKP2XotWhetjvghPqGDdwn0ULGRBw "Histogram & MA plot")

# Single gene Plot Counts
It can be useful to examine the counts of reads for a single gene across the groups. The function for making this plot is  _plotCounts_, which normalizes counts by sequencing depth and adds a pseudocount of 1/2 to allow for log scale plotting. 

The counts are grouped by the variables in  `intgroup`, where more than one variable can be specified. Here we specify the gene which had the smallest  _p_  value from the results table (res dataframe). You can select the gene to plot by rowname or by numeric index.

```r
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```

# PCA Plot

Visualise sample-sample distances with the Principal Components Analysis (PCA) 

 - Complementary approach to assess if samples have greater variance between experimental and control conditions than between replicates.
 - Aim is to **identify groups of features** (eg genes) that have something in common, such as expression patterns across different samples.
 - Result is principal components representing directions along which the variation in the originial multi-dimensional data is maximal, so that a few components (dimensions) can be used to represent thousands of mRNA data points.
 - Can visually represent variation of gene expression for different samples in a simple xy plot (instead of plotting thousands of genes per sample)/ Usually only the top 2 principal components (explaining the majority of the data variability) are displayed.
	 - identify unexpected patterns - batch effects; outliers
	 - does not identify unknown groupings

in R use `prcomp` function"

```r
#Load libraries
library(edgeR)
library(DESeq2)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)

pc = prcomp(t(rlog.norm.counts))
plot(pc$x[ ,1], pc$x[ ,2], col = colData(DESeq.ds)[ ,1], main = "PCA of seq.depth normlised\n and rlog-transformed read counts"
# PCA plot using DESeq2 based on ggplot2
P = plotPCA(DESeq.rlog) 
#plot cosmetics
P = P + theme_bw() + ggtitle("Rlog transformed counts") 
print(P)

data <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
```


![enter image description here](https://onlinecourses.science.psu.edu/stat857/sites/onlinecourses.science.psu.edu.stat857/files/lesson05/PCA_plot/index.gif)



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
# Mult-dimensional scaling
```r
#### Plot #8 - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries

#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together
#What pattern do we expect to see, given the types of libraries we have (technical replicates, biologal replicates, tumor/normal)?
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

# Calculate the differential expression results including significance
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
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
# Heatmap DE
```r
#### Plot #11 - Create a heatmap to vizualize expression differences between the eight samples
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes=results_genes[sig,"id"]
sig_gene_names=results_genes[sig,"gene_name"]
data=log2(as.matrix(gene_expression[sig_genes,data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,7), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names,col=rev(heat.colors(75)))
dev.off()

#The output file can be viewed in your browser at the following url:
#Note, you must replace __YOUR_IP_ADDRESS__ with your own amazon instance IP
#http://__YOUR_IP_ADDRESS__/workspace/rnaseq/de/ballgown/ref_only/Tutorial_Part3_Supplementary_R_output.pdf

#To exit R type:

quit(save="no")
```

# cummeRbund

Automatically generates many of the commonly used data visualisations including:
- distribution plots
- correlation plots
- MA plots
- Volcano plots
- Clustering, PCA & MDS plots (global relationships between conditions)
- Heatmaps
- gene/transcript level plots showing transcript structures & expression levels 

Input the output from Cufflinks, Cuffmerge & Cuffdiff

# Explore Read Count Patterns




Assess RNA-seq expression patterns with:
1. Pairwise Correlation
2. Hierarchical Clustering
3. Principal Component Analysis

# Pairwise correlation

- Pearson correlation coefficient, r is used to assess similarity between RNA-seq samples in a pair wise fashion.
- ENCODE guideline advises **>0.9** correlation should achieved for mRNA transcripts.
- in R use `cor( )` function

# Hierarchical clustering
- separate samples in an **unsupervised** fashion to see if samples of different conditions differ more than replicates within the same condition.
- pairwise compairsons of individual samples, grouped into "neighbourhoods" of similar samples.
- hierachical clustering analyses require decisions on:
	- how the disimilarity between pairs should be calculated?
	- how the disimilarity be used for clustering?
		- Pearson correlation coefficient
		- Euclidean distance (distance between two vectors) - calculate the distance using `linkage` function (complete, average, or single intercluster distance - complete intercluster distance is best and single IC distance is worst).

Can visualise DE genes as a Dendogram or Heatmap

## Clustered Heatmap
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

## Dendogram
Genes are sorted by adjusted p-value. Colours represent read counts.
- clustures obtained by cutting dendoram at a level where the jump between two nodes is large
- connected components form individual clusters
- clustering algorithms differ and there is no concensus on which is optimal
in R use `cor( )` `as.dist( )` and `hclust( )` to generate a dendogram




![enter image description here](http://www.sthda.com/english/sthda-upload/figures/cluster-analysis/009c-divisive-hierarchical-clustering-compute-diana-1.png)



# Explore normalised read counts
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

# Histogram

hist (DGE. results $ pvalue, col = " grey ", border = " white ", xlab = "", ylab = "", main = " frequencies of p- values ") # histogram of p-valus




# Variance Shrinkage

`DSeq2` and `edgeR` both offer means to reduce the variance using the dispersion mean trend using the entire dataset as a reference.

Low read counts that are highly variable will be assigned more homogenous read estimates --> variance resembles the majority of the genes and hopefully has a more stable variance

Regularise log-transformed values:
`DSeq.rlog = rlog(DESeq.ds, blind = TRUE)` #can set rlog to FALSE if there are large differences in a large proportion of the genes to avoid overestimating the dispersion
`rlog.norm.counts = assay(DESeq.rlog)`

`msd_plot = meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)` #show data on original scale and dont print plot
`msd_plot$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation)`


**SVD (singular value decomposition) analysis**

-   For doing this you can use the gene-level count table obtained from Kallisto. I wrote everything in R and I can send you some litterature which explains a bit the underlying math and idea. Also happy to speak about it over skype.

# Ballgown Visualisation

https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part2_ballgown.R
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEzOTI4NTkwNDEsMjEyMzc2NjIzMiwtND
Y0OTQ4NzE5LDk0Njg1MDg5MywtMzMwMjkwMTE5LDk1OTMyNzk4
OSwxODEwODI0NTQ2LC0xOTkwNjk3NDE1LDE0NDU0Nzk4MjMsOD
U5Njc3MjUzLDY4MDAxNjIxOCwxMzMwNjE1NzA4LDUzMDAxMDAw
NSwtODc2MDI1NTQ5LC0xMzk5NzM0NDA0LC0xMTE0NzY3NjIwXX
0=
-->