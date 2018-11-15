

> # RNA sequence protocol for assessing Alternative Splicing

# Contents

- This repository contains a protocol to analyse RNA-seq data, focusing on alternative splicing & polyadenylation, authored by Oliver Ziff. 
- The contents are based on multiple resources including:
	-  the [RNAseq worksheet](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf); 
	- [Biostars handbook](https://www.biostarhandbook.com/)
	- [Data Camp](https://www.datacamp.com/home)
	- [Coursera](https://www.coursera.org/specializations/bioinformatics)
	- and most importantly the experience of established experts in RNAseq analysis within [the Luscombe lab - my host laboratory](https://www.luscombelab.org/crickmembersdetail). 
- The protocol utilises a combination of bash `unix` commmand line and `R` scripts.

## Chapters
1. RNA seq workflow
2. Wet-lab RNA sequencing phase
3. Accessing sequencing data
4. QC of sequencing files
5. Alignment
6. Visualisation in IGV browser
7. QE of aligned reads
8. Read quantification
9. Differential expression analysis
10. Splicing analysis
11. Gene enrichment analysis 

# RNA-seq Workflow

## Wet-lab sequencing phase:
1. Extract & isolate RNA
2. Prepare library: break RNA into small fragments, convert to dsDNA, add sequencing adapters, PCR amplify
3. Strand Sequence the cDNA library: flow cell, base calling & quality score, replicates (technical = multiple lanes in flow cell; biological = multiple samples from each condition)
![Preparing RNA seq library](https://lh3.googleusercontent.com/RYpyReGfJbJOWjm20hzclqR6KUMkacZ6p_xaKvQs3piOTfxXdRiXUmiKAd45nHWj30cxJPVXmqTfnQ)
![enter image description here](https://lh3.googleusercontent.com/EBRN0O87F248JvjOzL_yHF1U328THjmXywtF4shxKxmzIwePgU-XR6ETv9Q0LCFP7bEcltsTXrN9hg)

## Bioinformatic phase:
https://www.biostarhandbook.com/rnaseq/rnaseq-intro.html

1. Process raw Reads: FATQ files download SRA, quality scores (Phred), paired vs single end sequence, FASTQC quality control, variability, spike-ins, blocking & randomise, filter out low quality reads & artifacts (adapter sequence reads).
2. Align (map) reads to reference genome (FASTA, GFF, GTF): annotation file (BED), alignment program (STAR, HISAT), reference genomes (GenCODE, Ensemble), generate genome index, create & manipulate BAM/SAM files containing sequence alignment data
3. Visualise & explore alignment data in IGV and R studio: ggplot2, bias identification QoRTs, 
4. Estimate Read Quantification (abundance) with gene based read counting 
5. Compare abundances between conditions & replicates (differential expression): Normalise, adjust each gene read counts for the total aligned reads  within each sample. Summarise data with pairwise correlation, hierarchical clustering, PCA analysis - look for differences between samples & identify outliers to consider excluding.

![Compare mutant vs wild type gene expression](https://lh3.googleusercontent.com/VtBLKXVhTx_hwbUNxN59byRcd2Ums76QpdRmtHYGUSo2wiwi5MkDEld8Eej6Bgsiqo25kJ4vxwtxNw)

![enter image description here](https://ycl6.gitbooks.io/rna-seq-data-analysis/Workflow.png) 
![enter image description here](https://www.rna-seqblog.com/wp-content/uploads/2016/02/typical.jpg)

## Requirements

On the CAMP cluscd ter most packages are preinstalled but to use them you need to use the module load function:
`ml STAR`
`ml ncbi-vdb`
`ml fastq-tools`
`ml SAMtools`
`ml RSeQC`
`ml QoRTs`
`ml multiqc`
`ml Subread`
`ml Java`
Use `module spider` to search for packages.

Install `conda` and activate `bioconda`

**Installing packages in R**
`install.package("package name")`
**Bioconductor** is a free software project for genomic analyses based on R programming. 
[Install Bioconductor](https://www.bioconductor.org/install/)
[Source]("https://bioconductor.org/biocLite.R")
`source ("https://bioconductor.org/biocLite.R")` 
`biocLite (“package_name“)`
`biocLite("erccdashboard")` # erccdashboard (for artificial spike in quantification) 
`biocLite("DESeq")`

Even though packages have been installed into R locally, then need to be brought into the working memory before using them:
`library("erccdashboard")`
`library("DESeq")`

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTA2NTQ4NDc4MCwtMTYxMDE0OTkyM119
-->