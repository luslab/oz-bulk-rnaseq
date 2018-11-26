


> # Visualise Differential Expression

# IGV

Firstly, visualise the most significantly DE genes in IGV
1. On local terminal `cd ~/bin/IGV_2.4.14/lib` & run IGV via command line on local terminal: `java -Xmx750m -jar igv.jar`
2. Set reference genome to Human (hg38) top left box.
3. Click File load from file > click Desktop > mount CAMP locally > click relevant BAM & BAI files (can load multiple at once).

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

Check that basic global patterns are met:
- do replicates have similar expression patterns; 
- do experimental conditions have differences
- For gene which you have prior knowledge about, you should check to see if they behaved as expected - a knockout gene should be very strongly down-regulated in the DGE analysis.
- Map the ORF identifiers from the read count matrix to the gene name --> retreive the rlog transformed read counts & log2 fold changes.

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

## Principal Components Analysis (PCA)

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

# MA plot
- provides global view of relationship between expression change in different conditions, average expression strength of genes and ability og algorithrm to detect differential expression
	- red dots are significant adjusted p<0.05
`plotMA(DGE.results, alpha = 0.05, main = "WT vs. SNF2 mutants", ylim = c(-4, 4))`

![enter image description here](https://lh3.googleusercontent.com/PdRsM9aHl3MTvEMKCYjYKQysVZ9MKxk943_XZ_JLLtAH0jTgZXKP2XotWhetjvghPqGDdwn0ULGRBw "Histogram & MA plot")

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
eyJoaXN0b3J5IjpbMTIwNzUzNDQwNCwtMTExNDc2NzYyMF19
-->