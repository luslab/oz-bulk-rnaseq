

> # RNA sequence protocol assessing for Alternative Splicing & Polyadenylation

- This repository contains a protocol to analyse RNA-seq data, focusing on alternative splicing & polyadenylation, authored by Oliver Ziff. 
- The contents are based on multiple resources including the [RNAseq worksheet](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf); [Biostars handbook](https://www.biostarhandbook.com/); [Data Camp](https://www.datacamp.com/home); [Coursera](https://www.coursera.org/specializations/bioinformatics); and most importantly the experience of established experts in RNAseq analysis within [my host laboratory](https://www.luscombelab.org/crickmembersdetail). 
- The protocol utilises a combination of bash `unix` commmand line and `R` scripts.

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
 











 














# Read Quantification
We first use RNA seq to determine the abundance of mRNA (cDNA) fragments, rather than the composition of the fragments. 

Different ways to quantify mRNA abundances:
1.  Counts: The number of reads overlapping with a transcript.
2.  RPKM/FPKM: Reads/Fragments per kilobase of transcript per millions of read mapped.
3.  TPM: Transcripts per million

## Gene-based read counting
- To compare the expression rates of individual genes between samples you need to **quantify the number of reads per gene.**
- Essentially you are **counting the number of overlapping reads**
- Need to clarify:
	- Overlap size (full read vs partial overlap)
	- Multi-mapping reads
	- Reads overlapping multiple genomic features of the same kind
	- Reads overlapping introns

## Spike-in control

- Spike-ins provide a known concentrations of transcripts that we can compare to the experimental samples
- A common spike in product is [ERCC ExFold RNA spike-in control mix](http://data.biostarhandbook.com/rnaseq/ERCC/ERCC-information.pdf) which is added to the experimental samples
- Fold change in transcript expression between 2 samples tells you about the difference between the 2; not about whether they are highly or lowly expressed.
- At lower transcript expression levels accuracy in determining fold change deteriorates. 

## Tools to count reads

- [htseq-count](http://htseq.readthedocs.io/en/release_0.10.0/index.html) has 3 modes union, intersection strict, and intersection nonempty (image above). 
- `featureCounts` counts reads if any overlap is found with a gene. Can exclude multi-overlap reads or include then for each gene that is overlapped. This is a package of Subread so need to `ml Subread` - Biostars advise this.
- `QoRTs` also does counting - Nobby uses this.

## featureCounts Workflow
ml Subread

Input:
- Gene feature file
- BAM file 

1. Count reads (estimate abundance) per sample:
```bash

# Create output folder
mkdir -p featureCounts

#set gene coordinates
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf
#set BAM input file
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/SRR5483788_Aligned.sortedByCoord.out.bam
#set Counts.txt output file
COUNTS=/home/camp/ziffo/working/oliver/projects/airals/featureCounts/D7_samples/counts_SRR5483788.txt

#run featureCounts command - by default it uses gene_id in the GTF file. Override this with gene_name attribute.
featureCounts -a $GTF -g gene_name -o counts.txt $COUNTS $BAM
```
This script produces 1 txt file per BAM file.
Using the * wildcard you can list all BAM files into 1 text file.

`featureCounts -a $GTF -g gene_name -o counts.txt /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/SRR5*_Aligned.sortedByCoord.out.bam`

The output file contains a column for each sample. 

For each BAM file there are 2 output files:
- featureCounts_results.txt has actual read counts per gene - tab delimited file where the first six columns contain feature specific information and the rest of the columns contain the read counts that overlap with that feature.
- featureCounts_results.txt.sumary gives quick overview of how many reads were assigned to genes. 

2. Find sequences with highest abundance

To find sequences with most hits sort by column 7: `cat counts.txt | sort -rn -k 7 | head`
Output table is in columns as:
```
Geneid            Chr         Start     End  Strand   Length  Hits
```

# Differential Gene Expression Analysis

**Tools for DGE:**
`edgeR` (better for false positives, less conservative, recommended if <12 replicates)
`DESeq`
`DESeq2` sample wise size factor
`limma-voom`
`cuffdiff`slow, cant support multifactored experiments, can detect differential isoforms, high false positives

## Compare replicates between conditions (i.e. differential expression)

Using R assess differential expression with DESeq2 using the output from featureCounts
```bash
ml R

#print out only columns representing Gene ID & sample abundances (ie remove Chr, Start, End, Strand & length)
cat counts.txt | cut -f 1,7-12 > simple_counts.txt

#download R script to run DESeq & DESeq2 using the simple_counts.txt file
curl -O http://data.biostarhandbook.com/rnaseq/code/deseq1.r
curl -O http://data.biostarhandbook.com/rnaseq/code/deseq2.r

#pass simple_counts.txt through the script specifying the design of the experiment 
## in this case = 3 x 3 (3 cases, 3 controls)
cat simple_counts.txt | Rscript deseq1.r 3x3 > results.txt
```
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
#Sort data by gene ID to paste into columns & select only foldchange and log2FoldChange. The results.txt file is already sorted according to padj
cat results.txt | sort | cut -f 1,5,6 > table

#How many genes are significantly differentially expressed (i.e. padj < 0.05)?
cat results.txt | awk ' $8 < 0.05 { print $0 }' > diffgenes.txt

#How many differentially expressed genes do we have?
cat diffgenes.txt | wc -l
```
Visualise the most significantly differentially expressed genes in IGV:
1. On local terminal `cd ~/bin/IGV_2.4.14/lib` & run IGV via command line on local terminal: `java -Xmx750m -jar igv.jar`
2. Set reference genome to Human (hg38) top left box.
3. Click File load from file > click Desktop > mount CAMP locally > click relevant BAM & BAI files (can load multiple at once).

## Gene Enrichment
https://www.biostarhandbook.com/ontology/gene-set-erichment.html
Now that you have identified the differentially expressed genes you need to identify their function.

**[Sequence Ontology](http://www.sequenceontology.org/browser/obob.cgi)**
There are >2,400 terms associated with sequences in the genome. Sequence Ontology defines sequence features used in biological annotations.
To search a defined sequence term use the [Sequence Ontology Browser](http://www.sequenceontology.org/browser/obob.cgi)

**[Gene Ontology](http://geneontology.org/)**
Connects each gene to one or more functions.
3 sub-ontologies for each gene product:
- Cellular Component (CC): cellular location where product exhibits its effect
- Molecular function (MF): How does gene work?
- Biological Process (BP): What is the gene product purpose?
Searching GO: use http://geneontology.org/ or https://www.ebi.ac.uk/QuickGO/

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

## Using pseudo-alignment to quantify transcript abundances

This is an alternative method to perform DE involving identifying which transcript that the reads originates from. Thus it aligns to the transcriptome (not to the genome).
This makes processing much quicker as it bypasses alignment to the genome but will not identify new transcripts.

Tools:
- [Kallisto](https://www.nature.com/articles/nbt.3519)
- Salmon

### Kallisto workflow
ml kallisto

```bash
#set shortcuts
REF=PATH_TO_FASTA_FILE.fa
IDX=PATH_TO_INDEX.idx
R1=PATH_TO_FASTQ_forward_strand.fq
R2=PATH_TO_FASTQ_reverse_strand.fq

#build kallisto index
kallisto index -i $IDX $REF

#create directory for kallisto data and set subdirectory for output called out
mkdir -p kallisto
OUTDIR=out
#run kallisto quantification with quant command. -o sets output directory. -b specifies the bootstap sample number.
kallisto quant -i $IDX -o $OUTDIR -b 100 $R1 $R2
```
You can run this for multiple samples as a For Loop (example [here](https://www.biostarhandbook.com/rnaseq/rnaseq-griffith-kallisto.html))
```bash
# Create output folder
mkdir -p kallisto

# Exit this script on any error.
set -euo pipefail

# This is the path & name of the reference transcriptome (not genome).
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.fa

# This the path & name of the index to build
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
        kallisto quant --single -l 200 -s 0.1 -i $IDX -o $OUTDIR -b 100 $R

        # Copy the abundance file to a proper name - i.e. remove the long path name so that it only contains information on the sample e.g. VCP_3.tsv
        cp $OUTDIR/abundance.tsv $OUTDIR.counts.tsv
    done
done

#concatenate all counts into 1 file that can be used for DESeq DE analysis.
paste /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/VCP_*.tsv  /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/CTRL_*.tsv | cut -f 1,4,9,14,19,24,29 > counts.txt
```
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


**Preparing an annotation:**
To **assess differential expression of exons**, create an annotation file where overlapping exons of different isoforms are split before running featureCounts. Use `dexseq_prepare_annotation.py` script of DEXSeq package or `QoRTs`.

# Gene Isoform counting

Gene isoforms are mRNA produced from the same locus but with diferent protein codeing sequences:

5 modes of alternative splicing are recognized:

1.  Exon skipping or cassette exon.
2.  Mutually exclusive exons.
3.  Alternative donor (5') site.
4.  Alternative acceptor (3') site.
5.  Intron retention.

![Alternative Splicing](https://en.wikipedia.org/wiki/Protein_isoform#/media/File:Alternative_splicing.jpg)

Strictly you should quantify reads that originate from transcripts (rather than genes as a whole).  Simple count-based approaches underperform when determining transcript level counts as they disregard reads that overlap with more than one gene. If the genomic feature becomes a transcript rather than a gene it keeps many reads that would have been discarded.

Programmes to quantify isoforms:
`Cufflinks`
`RSEM`
`eXpress`

These use a **deBruikin graph** to assign reads to an isoform if they are compatible with that transcript structure.
![enter image description here](https://www.frontiersin.org/files/Articles/169488/fgene-06-00361-r2/image_m/fgene-06-00361-g002.jpg)

Schema of a simple deBruijn graph-based transcript assembly. (A) Read sequences are split into (B) all subsequence k-mers (here: of length 5) from the reads. (C) A deBruijn graph is constructed using unique k-mers as the nodes and overlapping k-mers connected by edges (a k-mer shifted by one base overlaps another k-mer by k􀀀1 bases). (D) The transcripts are assembled by traversing the two paths in the graph

### Alternative approach to read counting
Ignore exactly where within a transcript a read originates from. Instead **focus on which transcript the read represents**. This approach does not generate a BAM file (alignment file) but instead produce a measure of how many reads indicate the presence of each transcript.

Tools for this approach:
`Sailfish` and more updated version `Salmon`
`Kallisto`

This approach is much faster than alignment-counting routines but **cant detect novel isoforms**. However, instead of direct isoform quantification, you can glean more accurate answers from alternative approaches, e.g., quantification of exons (Anders et al., 2012) or estimates of alternative splicing events such as exon skipping, intron retention etc. (e.g., MISO [Katz et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037023/), rMATS [Shen et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280593/)).

The main limitations to assigning reads to transcripts are:
- annotation transcripts are inconsistent
- many isoforms with very different lengths
- anti-sense and overlapping transcripts of different genes


# Differential Gene Expression Analysis

1.  _Within-sample_  comparisons: compare the expression of genes within the same experiment eg in this experiment does gene A express at a higher level than gene B?
2.  _Between-sample_  comparisons: compare the expression of a gene between experimental conditions aka pairwise comparison e.g. has the gene expression for gene A  changed across different experimental conditions? if yes, then differentially expressed (DE).

## Normalising and Log Transforming Read Counts
The number of sequenced reads depends on:
1. expression level
2. read length
3. sequencing depth
4. expression of all other genes in the sample

To compare two conditions then you must look at the fraction of transcripts assigned to a specific gene over the total number of reads which differs drastically between samples.

Normalisation is performed to ensure that systematic effects not related to the biological differences between samples are removed. 

Methods of normalising:
- total count
- Counts/million
- **DSeq size factor** using R
- TMM (trimmed mean of M values)
- upper quartile


**`DSeq`:**

 - Takes the `featureCounts` (raw read counts) --> read them into R --> normalise for sequencing depth differences
 
`DSeqDataSet` contains all the information in R. 
`rowRanges( )` rows = variables (genes, transcripts, exons)  - has info about genes (chr, start, end, strand, ID)
`rownames` is the unique sample names
`colData( )` columns = samples 
`assay( )` stores count data with genes and samples. similar to `countData`

**DSeq process:**
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
8. Generate the DSeqDataSet
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


### Log Transformation of Sequencing Depth Normalised read counts

Most downstream analyses work best on log scales of read counts. Usually *log2* but occasionally *log10*.

Log2 transform read counts:
`log.norm.counts = log2(counts.sf_normalized + 1)` #use a pseudocount of 1

### Plot images to visually explore normalised read counts:
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

### Variance Shrinkage

- `DSeq2` and `edgeR` both offer means to reduce the variance using the dispersion mean trend using the entire dataset as a reference.
- Low read counts that are highly variable will be assigned more homogenous read estimates --> variance resembles the majority of the genes and hopefully has a more stable variance

Regularise log-transformed values:
`DSeq.rlog = rlog(DESeq.ds, blind = TRUE)` #can set rlog to FALSE if there are large differences in a large proportion of the genes to avoid overestimating the dispersion
`rlog.norm.counts = assay(DESeq.rlog)`

`msd_plot = meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)` #show data on original scale and dont print plot
`msd_plot$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation)`

### Explore Global Read Count Patterns
- check that basic global patterns are met:
	- do replicates have similar expression patterns; 
	- do experimental conditions have differences

3 most commonly ways to **assess RNA-seq expression patterns**:

**1. Pairwise correlation**

- Pearson correlation coefficient, r is used to assess similarity between RNA-seq samples in a pair wise fashion.
- ENCODE guideline advises **>0.9** correlation should achieved for mRNA transcripts.
- in R use `cor( )` function

**2. Hierarchical clustering**
- separate samples in an **unsupervised** fashion to see if samples of different conditions differ more than replicates within the same condition.
- pairwise compairsons of individual samples, grouped into "neighbourhoods" of similar samples.
- hierachical clustering analyses require decisions on:
	- how the disimilarity between pairs should be calculated?
	- how the disimilarity be used for clustering?
		- Pearson correlation coefficient
		- Euclidean distance (distance between two vectors) - calculate the distance using `linkage` function (complete, average, or single intercluster distance - complete intercluster distance is best and single IC distance is worst).
- Output = dendogram:
	- clustures obtained by cutting dendoram at a level where the jump between two nodes is large
	- connected components form individual clusters
	- clustering algorithms differ and there is no concensus on which is optimal
in R use `cor( )` `as.dist( )` and `hclust( )` to generate a dendogram
code on page 55

![enter image description here](http://www.sthda.com/english/sthda-upload/figures/cluster-analysis/009c-divisive-hierarchical-clustering-compute-diana-1.png)

**3. Principal Components Analysis (PCA)**

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

# Differential Gene Expression (DGE) Analysis
the 2 tasks of DGE are to:
1. Estimate the magnitude **fold change** of differential expression in read counts between different conditions, accounting for differences in sequencing depth & variability
2. Estimate the **significance** of the difference, accounting for multiple testing

Null hypothesis = the mean read counts of genes in different samples are equal between different conditions.

Model read counts using *Poisson distribution*. This is useful as:
- individual reads can be interpreted as binary
- we can model the discrete probability distribution of the number of reads identified in the sequenced library
- the pool of possible reads is huge, but the proportion of reads belonging to gene x is small
- variance = mean in Poisson distribution --> from the mean read count per condition, we know the variance --> we can identify genes with greater differences between conditions
- repeat libraries can be well approximated using poisson - biological replicates have relatively high variance and this *over-dispersion* can be captured with the *negative binomial* distribution

Need to estimate the **mean** and **dispersion** from the read counts:
- precision depends on the number & variation of replicates (i.e. the statistical power). 
- For RNA-seq typically there is only 2-3 replicates creating poor precision. There are tools to compensate for this by looking across genes with similar expression to reduce a given genes variance towards the regressed values.


![enter image description here](https://lh3.googleusercontent.com/LVvCl3GXhNzUx5lyTrHsr0z_ZmI0nb51TBiY1-53VifMuYW8HR9-X54sfLwoH5gFyqahHOm8_QaWhg "Comparison of DGE programs")

##  DE analysis in R

1. Normalise matrix between genes 
2. Plot the normalised matrix

Clustered heatmap 

**DSeq2 workflow**
This is *performed on the raw read* counts and not the transformed normalised reads.
`str(colData(DESeq.ds)$condition)` #use the levels of condition to determine the comparison order
`colData(DESeq.ds)$condition = relevel(colData(DESeq.ds)$condition, "WT")` #set wild type as the first level factor (the mutants should be compared with the control - use the wildtype as the denominator)
`DESeq.ds = DESeq(DESeq.ds)` #run the DGE analysis
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

### Exploratory Plots

**Histograms**
hist (DGE. results $ pvalue, col = " grey ", border = " white ", xlab = "", ylab = "", main = " frequencies of p- values ") # histogram of p-valus

**MA plot**
- provides global view of relationship between expression change in different conditions, average expression strength of genes and ability og algorithrm to detect differential expression
	- red dots are significant adjusted p<0.05
`plotMA(DGE.results, alpha = 0.05, main = "WT vs. SNF2 mutants", ylim = c(-4, 4))`

![enter image description here](https://lh3.googleusercontent.com/PdRsM9aHl3MTvEMKCYjYKQysVZ9MKxk943_XZ_JLLtAH0jTgZXKP2XotWhetjvghPqGDdwn0ULGRBw "Histogram & MA plot")

**Heatmaps**
- show expression values across individual samples
many R functions to do this:
`aheartmap( )`
`gplots::heatmap.2( )`
`pheatmap::pheatmap( )`
page 60 for R code

![enter image description here](http://bioinfo.cipf.es/babelomicstutorial/_media/images:differential_expression_example:heatmap.png)
Genes are sorted by adjusted p-value. Colours represent read counts.

**Read counts of single genes**
- For gene which you have prior knowledge about, you should check to see if they behaved as expected. For example a knockout gene should be very strongly downregulated in the DGE analysis.
- Map the ORF identifiers from the read count matrix to the gene name --> retreive the rlog transformed read counts & log2 fold changes.
- Use an annotation database within R specific for your sample eg for human use org.Hs.eg.db https://www.bioconductor.org/packages/release/data/annotation/
p61 for R code
Edge R workflow: page 62





@Raphaelle used: 
 **Splicing analysis**  
   
-   I first used  [VAST-tools](https://github.com/vastgroup/vast-tools) which performs alignment for you. So basically you submit your fastq files directly. Have a look at the GitHub vignette as it is rather complete however please do not hesitate to contact me if you want help with shell scripting as you will need to run this as a loop. Or Nobby will certainly be happy to help on CAMP (I am working from UCL cluster and have never logged onto CAMP).
-   Then to perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.

  **3' UTR isoforms analysis**
-   I have written an entire pipeline for this which I can explain and share with you scripts when needed. But basically I extract genome-wide coverage using bedtools, then extract regions of continuous coverage along genome, then intersect these with Ensembl annotated regions, extend 3' UTR. Finally to annotate all alternative 3' UTR isoforms I then run an algo which identifies shifts in coverage along 3' UTR which are expected to occur at PAS sites.

 **Differential gene expression**

-   I first use  [Kallisto](https://pachterlab.github.io/kallisto/)  which is a really user-friendly algo which extract both gene and transcript level gene expression directly from fastq files. So here again I directly used the raw fastq files.
-   The I used  [Sleuth](https://pachterlab.github.io/sleuth/)  (also developed by Pachter lab) to perform differential gene and transcript expression analysis.

**SVD (singular value decomposition) analysis**

-   For doing this you can use the gene-level count table obtained from Kallisto. I wrote everything in R and I can send you some litterature which explains a bit the underlying math and idea. Also happy to speak about it over skype.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMjEwMjY2MjE3MiwtMTY4ODMwNjA2NF19
-->