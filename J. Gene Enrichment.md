> # Gene Enrichment

https://www.biostarhandbook.com/ontology/gene-set-erichment.html
Now that you have identified the differentially expressed genes (either up or down regulated), we need to identify their function. This is performed in a Gene Set Overlap analysis. Do the up or down regulated genes have anything in common?  The hypergeometric test (Fisher's test) is used to test whether genes enriched in the list of up or down regulated genes occur more often than expected by chance.

Gene Ontology (GO) is a vocabulary used to describe genes functions with respect to 3 aspects:
- Cellular Component (CC): cellular location where product exhibits its effect
- Molecular function (MF): How does gene work?
- Biological Process (BP): What is the gene product purpose?

Searching GO: use http://geneontology.org/ or https://www.ebi.ac.uk/QuickGO/


# GO analysis
https://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

Use reporting tools to write a table of GO analysis results to a HTML file. 
Select genes of interest > run hyperGTest > make GO report

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
GenTable( myGOdata, goTestResults )
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
```

# Plot GO p-values as Barplot
https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

Always filter out redundant term and select the 10 most significant otherwise with barplot you can't read anything
I guess the purpose of the above package GOplot is to plot many more terms than 10.
```r
#If you want only to plot the GO terms in one condition
dat        <- enrich$value # the -log10(P-value)
names(dat) <- enrich$Terms #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1)

#If you want to plot side by side the same GO term but for 2 different conditions
dat               <- -cbind(enrich.bp$p.value.norm.ngf,enrich.bp$p.value.norm.nt3)[match(least.transported,enrich.bp$names.goterms.),]
rownames(dat)    <- enrich.bp$Term[match(least.transported,enrich.bp$names.goterms.)]
dat              <- dat[order(dat[,1],-dat[,2]),]
mycols           <- c("#81A4D6","#AF71AF")
L                <- nrow(dat)
val1             <- -as.vector(dat[,1])
val2             <-  as.vector(dat[,2])
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
plot(c(L+3,0),xlim=c(min(val1)-10,max(val2)),type = "n",frame=F,yaxt="n",ylab="")
mp=barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col=mycols[2])
barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col=mycols[1])
text(x=(min(val1)-5),y=mp,lab=rownames(dat),cex=0.6)
mtext(side=1,line=2,text="Z-score",cex=0.6)
abline(v=c(-2.5,2.5),col="grey",lty=2)
```

```r
#If you want only to plot the GO terms in one condition
dat        <- enrich$value # the -log10(P-value)
names(dat) <- enrich$Terms #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1)



#If you want to plot side by side the same GO term but for 2 different conditions
dat               <- -cbind(enrich.bp$p.value.norm.ngf,enrich.bp$p.value.norm.nt3)[match(least.transported,enrich.bp$names.goterms.),]
rownames(dat)    <- enrich.bp$Term[match(least.transported,enrich.bp$names.goterms.)]
dat              <- dat[order(dat[,1],-dat[,2]),]

mycols           <- c("#81A4D6","#AF71AF")
L                <- nrow(dat)
val1             <- -as.vector(dat[,1])
val2             <-  as.vector(dat[,2])
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
plot(c(L+3,0),xlim=c(min(val1)-10,max(val2)),type = "n",frame=F,yaxt="n",ylab="")
mp=barplot(height = val2,add = TRUE,axes = FALSE,horiz=T,col=mycols[2])
barplot(height = val1,add = TRUE,axes = FALSE,horiz=T,col=mycols[1])
text(x=(min(val1)-5),y=mp,lab=rownames(dat),cex=0.6)
mtext(side=1,line=2,text="Z-score",cex=0.6)
abline(v=c(-2.5,2.5),col="grey",lty=2)
```

# BLAST
Basic Local Alignment Search Tool

Use blastx for RNA

[Swiss-Prot database](https://www.uniprot.org/):
-includes TrEMBL. Swiss-prot is a manually curated subset of this. 
Gives information on Function & GO on Proteins & Genes.

# **[Sequence Ontology](http://www.sequenceontology.org/browser/obob.cgi)**
There are >2,400 terms associated with sequences in the genome. Sequence Ontology defines sequence features used in biological annotations.
To search a defined sequence term use the [Sequence Ontology Browser](http://www.sequenceontology.org/browser/obob.cgi)

**[Gene Ontology](http://geneontology.org/)**



GO Download page: http://geneontology.org/page/download-annotations

2 files to download: 
1. definition (term) file `wget http://purl.obolibrary.org/obo/go.obo`
2. association file `wget http://geneontology.org/gene-associations/goa_human.gaf.gz`. In GAF compressed format defined at http://geneontology.org/page/go-annotation-file-gaf-format-21
Contains both gene and protein IDs.

To search the function of a gene use the [GeneCards](http://www.genecards.org/) database to easily locate the gene by name.

### Gene Set Enrichment Analysis
- Identify common characteristics within a list of genes. When using GO terms, this is called "functional enrichment"
- Most common variant is the ORA (over-representation analysis): 
	- examines genes in a list > 
	- summarises the GO annotations for each gene > 
	- determines if any annotations are statistically over-represented.

2 x 2 table:

|  | Differentially Expressed | Not Differentially Expressed
|--|--|--|
| Gene Ontology Term | 50 | 200
| No Gene Ontology Term| 1950 | 17800


 **GO enrichment tools** 
 The best are:
- [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) Bioconductor package (used by Luisier et al 2018)
- AgriGO Web-based GO Analysis Toolkit and Database for Agricultural Community.
- DAVID This is the GO tool biologists love. It is the "most generous" of them all, as it produces copious amounts of output. 
- Panther is offered directly from the GO website. It produces very limited information and no visualization.
- goatools a command line Python scripts.
- ermineJ Standalone tool with easy to use interface. It has detailed documentation.
- GOrilla Web-based; good visualization; downloadable results (as Excel files); easily see which genes contribute to which enriched terms; results pages indicate date of last update GO database (often within the last week).
- ToppFun - a suite of tools for human genomes that integrates a surprising array of data sources.

#### Gene Set Enrichment Analysis Methodology
1. Enter gene names to be studies
2. Enter background gene names (usually all genes for the organism)
3. Perform statistical comparison

[g:Profiler](https://biit.cs.ut.ee/gprofiler/) performs functional enrichment analysis and analyses gene lists for enriched features.  Very good visualiser of GO terms. 
[g:sorter](https://biit.cs.ut.ee/gprofiler/gsorter.cgi) finds similar genes in public  transcroptomic data. Input = single gene & dataset of interest. Result = sorted gene list similarly expressed with gene of interest. For global gene expression analyses, across different species use [Multi Experiment Matrix](https://biit.cs.ut.ee/mem/) tool.

The **Functional Annotation Tool** maps the genes to annotation content providing a summary of the biological interpretation of the data.

Perform **Fisher's Exact Test** to measure gene enrichment in annotation terms. The EASE score is a slightly modified Fisher's Exact p-value. The smaller to p-value, the more enriched the term.

# Trinity Platform

https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis

Trinotate web
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTk1MzQwNzA5NiwxMjY5ODI3MDE4LDkyNj
MyOTE0MSwtMTk0NTc1ODY4NSw4NDIzNzA4MzQsMTg1NjE0MjE3
MSwtMTE1MjQwNjMsMTc0NDQ3NjUxOCwtMTcxMzQ4MjI2OCwxMz
IxMjE1OTA3LDk0NzUyMDQ4OCw3OTc5NDUwMTcsNDg4NDU3Nzc3
LC05NDIwMTQzMCwxNTI4NTgxNTkzXX0=
-->