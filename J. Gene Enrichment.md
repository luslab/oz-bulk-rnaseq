> # Gene Enrichment

https://www.biostarhandbook.com/ontology/gene-set-erichment.html
Now that you have identified the differentially expressed genes (either up or down regulated), we need to identify their function. This is performed in a Gene Set Overlap analysis. Do the up or down regulated genes have anything in common?  The hypergeometric test (Fisher's test) is used to test whether genes enriched in the list of up or down regulated genes occur more often than expected by chance.

Gene Ontology (GO) is a vocabulary used to describe genes functions with respect to 3 aspects:
- Cellular Component (CC): cellular location where product exhibits its effect
- Molecular function (MF): How does gene work?
- Biological Process (BP): What is the gene product purpose?

Searching GO: use http://geneontology.org/ or https://www.ebi.ac.uk/QuickGO/


# GO analysis

Use reporting tools to write a table of GO analysis results to a HTML file. 
Select genes of interest > run hyperGTest > make GO report
```r
library(GOstats)
library(topGO)
library(org.Mm.eg.db)

# subset results table to only genes with sufficient read coverage
resTested <- resLFC1[ !is.na(resLFC1$padj), ]
# construct binary variable for UP or DOWN regulated
genelistUp <- factor( as.integer( resTested$padj < .1 & resTested$log2FoldChange > 0 ) )
names(genelistUp) <- rownames(resTested)

### Test UPREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistUp, nodeSize = 10,
   annot = annFUN.org, mapping = "org.Hs.eg.db", ID="ensembl" )
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


### Test DOWNREGULATED GENES





tt <- topTags(edgeR.de, n = 1000, adjust.method = 'BH', sort.by = 'p.value')
selectedIDs <- rownames(tt$table)
universeIDs <- rownames(mockRnaSeqData)
goParams <- new("GOHyperGParams", + geneIds = selectedIDs, + universeGeneIds = universeIDs, + annotation ="org.Mm.eg" , + ontology = "MF", + pvalueCutoff = 0.01, + conditional = TRUE, + testDirection = "over")
goResults <- hyperGTest(goParams)

goReport <- HTMLReport(shortName = 'go_analysis_rnaseq', + title = "GO analysis of mockRnaSeqData", + reportDirectory = "./reports")
publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Mm.eg", + pvalueCutoff= 0.05)
finish(goReport)
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
eyJoaXN0b3J5IjpbMTc0NDQ3NjUxOCwtMTcxMzQ4MjI2OCwxMz
IxMjE1OTA3LDk0NzUyMDQ4OCw3OTc5NDUwMTcsNDg4NDU3Nzc3
LC05NDIwMTQzMCwxNTI4NTgxNTkzXX0=
-->