> # Gene Enrichment

https://www.biostarhandbook.com/ontology/gene-set-erichment.html
Now that you have identified the differentially expressed genes (either up or down regulated), we need to identify their function. This is performed in a Gene Set Overlap analysis. Do the up or down regulated genes have anything in common?  

# Gene Ontology
https://www.biostarhandbook.com/ontology/gene-ontology.html

Gene Ontology (GO) is a vocabulary used to describe genes functions with respect to 3 aspects:
- Cellular Component (CC): cellular location **where** product exhibits its effect
- Molecular function (MF): **How** does gene work?
- Biological Process (BP): **What** is the gene product purpose?

GO terms are orgnised in directed acyclic graphs where relationships and heirachy are represented. View GO relationships at http://geneontology.org/ or https://www.ebi.ac.uk/QuickGO/


# GO functional analysis
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002375
https://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

There is 2 steps to this:
### 1. Overrepresentation analysis ORA (using Fisher test)

Look at collection of genes within a pathway found in a gene set (rather than individual genes) 
Is the more differentially expressed genes within a gene collection (e.g. RNA metabolism) than we would expect if genes were expressed randomly.

### 2. Enrichment sets

Use reporting tools to write a table of GO analysis results to a HTML file. 
Select genes of interest > run hyperGTest > make GO report

The hypergeometric test (Fisher's test) is used to test whether genes enriched in the list of up or down regulated genes occur more often than expected by chance.
Fisher Exact test = is there a relationship between 2 categorical variables? 
2x2 contingency table: present in gene set vs member of gene set:

|  | Differentially Expressed | Not Differentially Expressed
|--|--|--|
| Gene Ontology Term | 50 | 200
| No Gene Ontology Term| 1950 | 17800

GO Download page: http://geneontology.org/page/download-annotations
2 files to download: 
1. definition (term) file `wget http://purl.obolibrary.org/obo/go.obo`
2. association file `wget http://geneontology.org/gene-associations/goa_human.gaf.gz`. In GAF compressed format defined at http://geneontology.org/page/go-annotation-file-gaf-format-21
Contains both gene and protein IDs.
To search the function of a gene use the [GeneCards](http://www.genecards.org/) database to easily locate the gene by name.

### Gene Set Enrichment Analysis
- Identify common characteristics within a list of genes. When using GO terms, this is called "functional enrichment"
- Most common variant is the ORA (over-representation analysis): 

Examine genes in a list > summarises the GO annotations for each gene > determines if annotations are statistically over-represented.

 **GO enrichment tools** 

- [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) Bioconductor package (used by Luisier et al 2018)
- AgriGO Web-based GO Analysis Toolkit and Database for Agricultural Community.
- DAVID This is the GO tool biologists love - produces copious amounts of output. 
- Panther is offered directly from the GO website. It produces very limited information and no visualization.
- goatools a command line Python scripts.
- ermineJ Standalone tool with easy to use interface. It has detailed documentation.
- GOrilla Web-based; good visualization; downloadable results (as Excel files); easily see which genes contribute to which enriched terms; results pages indicate date of last update GO database (often within the last week).
- ToppFun - a suite of tools for human genomes that integrates a surprising array of data sources.
- [g:Profiler](https://biit.cs.ut.ee/gprofiler/) performs functional enrichment analysis and analyses gene lists for enriched features.  Very good visualiser of GO terms. 
- [g:sorter](https://biit.cs.ut.ee/gprofiler/gsorter.cgi) finds similar genes in public  transcroptomic data. Input = single gene & dataset of interest. Result = sorted gene list similarly expressed with gene of interest. For global gene expression analyses, across different species use [Multi Experiment Matrix](https://biit.cs.ut.ee/mem/) tool.

# Workflow

#### Gene Set Enrichment Analysis Methodology
1. Enter gene names to be studies
2. Enter background gene names (usually all genes for the organism)
3. Perform statistical comparison

The **Functional Annotation Tool** maps the genes to annotation content providing a summary of the biological interpretation of the data.

Perform **Fisher's Exact Test** to measure gene enrichment in annotation terms. The EASE score is a slightly modified Fisher's Exact p-value. The smaller to p-value, the more enriched the term.

## DAVID
https://www.biostarhandbook.com/ontology/david-analysis.html
https://david.ncifcrf.gov/home.jsp

Step 1: Load a gene list into the [site](https://david.ncifcrf.gov/tools.jsp)
Using the output of DESeq2 resOrderedDT datatable can extract the Entrez terms.
Paste in only the significant events (padj < 0.05) or top 100 events

Step 2: Select Gene identifier
Entrez (or Ensembl depending on which column from the DT you extract)

Step 3: Clarify List Type
Specify Gene List as you only pasted in significant events (Background is if you paste in entire gene list)

Submit > Select input species (Homo sapiens)

Explore details through the Functional Annotation Tools:  table, chart & clustering reports
Export & save results



```r

### GO analysis
library(topGO)
library(GOstats)
library(goseq)
library(org.Mm.eg.db)
library(GO-lite)
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
```

## Filter out redundant terms 

Take resTestedDT datatable
Filter out the non-significant say padj<0.05
if there are more than 150 then select only the top 100.
Extract columns GO term & padj
Paste into Revigo box via Excel
To filter out redundant terms first run [Revigo](http://revigo.irb.hr/). 


Manually curate the list by looking at genes content in the GO as even with Revigo they can be very redundant. In case you don't know how to get the list, please let me know and I will send you an R little command.

```r
require("GO.db")
require("topGO")
require("biomaRt")
require("org.Hs.eg.db")

#How to get GeneSymbols
x                  <- org.Hs.eg.db
mapped_genes       <- mappedkeys(x)
GO2geneID          <- as.list(x[mapped_genes])
geneID2GO          <- as.list(org.Hs.eg.db)
goterms            <- Term(GOTERM)
geneNames          <- names(GO2geneID)
mart              <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'ensembl.org')
anno              <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name","entrezgene"), mart = mart)
myGS              <- anno$external_gene_name
eid               <- anno$entrezgene
eid               <- eid[!is.na(eid)]
geneList          <- factor(as.numeric(geneNames%in%as.character(eid)))
names(geneList)   <- geneNames
mysampleGO        <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = as.character(eid),nodeSize = 10,annot = annFUN.GO2genes,GO2gene=geneID2GO)


goID              <-"GO:0008088"
go.entrez         <- genesInTerm(mysampleGO, goID)
go.entrez         <- go.entrez[[1]]#To extract all genes related to this term
go.gs             <- unique(as.character(anno$external_gene_name)[which(anno$entrezgene%in%go.entrez)])
```

### Test UPREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_BP_table <- GenTable( myGOdata, goTestResults )
enrich.ur.bp <- transform(UR_BP_table,result1 = as.numeric(result1))
enrich.ur.bp$value <- -log10(enrich$result1) # add column with -log10 of P-value
### Plot GO p-values as as bar plot
dat        <- enrich.ur.bp$value # the -log10(P-value)
names(dat) <- enrich.ur.bp$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.ur.cc <- transform(UR_CC_table,result1 = as.numeric(result1))
enrich.ur.cc$value <- -log10(enrich$result1)
### Plot GO p-values as as bar plot
dat        <- enrich.ur.cc$value # the -log10(P-value)
names(dat) <- enrich.ur.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_MF_table <- GenTable( myGOdata, goTestResults )
enrich.ur.mf <- transform(UR_MF_table,result1 = as.numeric(result1))
enrich.ur.mf$value <- -log10(enrich$result1)
### Plot GO p-values as as bar plot
dat        <- enrich.ur.mf$value # the -log10(P-value)
names(dat) <- enrich.ur.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)


### TEST DOWNREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_BP_table <- GenTable( myGOdata, goTestResults )
enrich.dr.bp <- transform(DR_BP_table,result1 = as.numeric(result1))
enrich.dr.bp$value <- -log10(enrich$result1)
### Plot GO p-values as as bar plot
dat        <- enrich.dr.bp$value # the -log10(P-value)
names(dat) <- enrich.dr.bp$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.dr.cc <- transform(DR_CC_table,result1 = as.numeric(result1))
enrich.dr.cc$value <- -log10(enrich$result1)

### Plot GO p-values as as bar plot
dat        <- enrich.dr.cc$value # the -log10(P-value)
names(dat) <- enrich.dr.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_MF_table <-GenTable( myGOdata, goTestResults )
enrich.dr.mf <- transform(DR_MF_table,result1 = as.numeric(result1))
enrich.dr.mf$value <- -log10(enrich$result1)

### Plot GO p-values as as bar plot
dat        <- enrich.dr.mf$value # the -log10(P-value)
names(dat) <- enrich.dr.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)
```

# Plot GO p-values as Barplot
https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

Always filter out redundant term and select the 10 most significant otherwise with barplot you can't read anything
I guess the purpose of the above package GOplot is to plot many more terms than 10.
```r
### Plot GO p-values as as bar plot
# add column with -log10 of P-value
sapply(UR_BP_table, class)
enrich <- transform(UR_BP_table,result1 = as.numeric(result1))
enrich$value <- -log10(enrich$result1)

#If you want only to plot the GO terms in one condition
dat        <- enrich$value # the -log10(P-value)
names(dat) <- enrich$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,10,2,3),cex=0.7)
barplot(height = dat,horiz=T,las=1, font.size = 20)


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


<!--stackedit_data:
eyJoaXN0b3J5IjpbMTM4MTUyNjc3NSwxNDA3ODAxMTg4LDc1MT
MyODA0OCwxMDg4MTUyMDc1LC0xMjI5MzU4NzY4LDIwOTkxOTc4
MjYsMTI1MTk3MjYxOSwtMjAwNTI4ODc5NywtODg4Mjg2MTE4LD
EyMjU1NjE1NDgsMTExMzU4MDE2MiwtMTk2MDUyMjkwNywyMDY3
NTM3MDI0LDEzOTQxNTU0MzEsLTUxMjI1NDA1MiwtMTcxMDU1Nz
UwNSwxOTUzNDA3MDk2LDEyNjk4MjcwMTgsOTI2MzI5MTQxLC0x
OTQ1NzU4Njg1XX0=
-->