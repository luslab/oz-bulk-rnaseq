> # Quality Control of Aligned Reads
 The main STAR output file for QC analysis is **Aligned.sortedByCoord.out.bam**.

After aligning and before performing downstream analyses check for:
1. Excessive amounts of reads not aligned
2. Obvious biases in the read distributions: 3' & 5' bias, base distribution, nucleotide content
3. Read quality
4. Sequencing depth
5. Similarity between replicate samples

Typical Biases of RNA-seq alignment:
- many reads aligned to introns indicates: 
	- incomplete poly(A) enrichment 
	- abundant presence of immature transcripts
- many reads align outside of annoted gene sequences (intergenic reads) indicates:
	- genomic DNA contamination
	- abundant non-coding transcripts
- over representation of 3' portions of transcripts indicates RNA degradation

### Tools
1. RSeQC package http://rseqc.sourceforge.net/
2. [QoRTs package](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf)
3. MultiQC

# RSeQC
ml RSeQC
ml QoRTs
ml MultiQC
ml R

## Assess gene body coverage for 3' & 5' Bias

Use **[RSeQC](http://rseqc.sourceforge.net/)** `geneBody_coverage.py` script:
- First it divides each transcript into 100 sections (irrespective of the transcript length)
![enter image description here](http://rseqc.sourceforge.net/_images/geneBody_workflow.png)
- then counts reads overlapping with each section
- produces 2 plots showing abundance of reads (coverage) across transcript bodies. Ideally you want to see equal coverage of sections across the whole transcript length.
![enter image description here](http://rseqc.sourceforge.net/_images/Aug_26.geneBodyCoverage.curves.png)
- in this plot there are 2 groups of transcripts. 1 group are symmetrical and have good transcript coverage. The other group has a skew towards the 3' end. This skew comes from:
	- polyA reads
	- degradation of RNA > check the RIN numbers.

Assess each BAM file separately:
```bash
ml RSeQC
ml R

BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC

#run each BAM file into geneBody coverage package using For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="geneBody_coverage.py -i $SAMPLE -r $BED -o ${OUT}/${SRRID} -f pdf"
done
```

This ouptut is 2 figures (line graph & heatmap) to visualise 3' or 5' bias. Each BAM file is represented by a different line. If you detect 3' bias at this stage you can either re-sequence (costly) or adjust for this bias in downstream analysis.

## Estimate RIN/TIN
Evaluate RNA integrity at transcript level (transcript integrity number) - analogous to to RIN. Measure of RNA quality. 

![enter image description here](https://www.researchgate.net/profile/Benjamin_Sigurgeirsson/publication/260841079/figure/fig5/AS:296675668185106@1447744400111/Gene-body-coverage-on-average-for-each-group-Both-RIN-10-and-RiboMinus-show-even.png)
Colours represent different RIN values (RIN 0 = degraded; RIN 9 = high quality). The RIN 0 line (degraded RNA) shows more 3' bias.

Determine a measure of mRNA degradation in silico using RSeQCs tin.py script to produce a TIN.
- TIN 0 (worst) - 100 (best). TIN 60 = 60% of transcript has been covered.
```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed

#run each BAM file into tin.py using a For Loop
for SAMPLE in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="tin.py -i $SAMPLE -r $BED"
done
```
Output is an xls file and a summary txt file (mean; median; SD values across all genes in sample).
Visualise TIN as boxplots in [Rstudio](https://github.com/friedue/course_RNA-seq2015/blob/master/03_mRIN.R) using ggplot:
```R
## Compare the distribution of mRIN (= TIN) values across different samples

# Check the local working directory in Rstudio
getwd()
#move the tin.xls files from CAMP to this local working directory. 
#run the Rscript on the xls files from https://github.com/friedue/course_RNA-seq2015/blob/master/03_mRIN.R
#This script will result in boxplots based on mRIN values calculated by## RSeQC

# list the respective files
TIN.files <- list.files(pattern = "xls") 
# lapply will iterate over the list of file names in TIN.files, read the respective tables and return a list of data frames where each data frame corresponds to one of the original xls tables
TIN.list <- lapply(TIN.files, function(x) read.table(x, header = TRUE)[c("geneID", "TIN")]) 
# to give meaningful names to each data frame, we can use regex on the file names the regex will have to be adjusted if you use different files
names(TIN.list) <- gsub("UHR-(.*)_Aligned.*", "\\1", TIN.files) 
# make a long data frame that is suitable for ggplot2 plotting
TIN.df <- as.data.frame(do.call(rbind, TIN.list)) 
# add a column that indicates the sample type for each gene ID and TIN value,here, I use the information from the original data frame's name in the TIN.list which is kept in the row.names of TIN.df
TIN.df$sample <- gsub("\\.[0-9]+", "", row.names(TIN.df)) 
# make the boxplots
library(ggplot2)
ggplot(data = TIN.df, aes(x = sample, y = TIN)) + geom_boxplot(notch=TRUE) 
# excluding genes with TIN = 0, which are most likely due to lack of read coverage
ggplot(data = subset(TIN.df, TIN > 0), aes(x = sample, y = TIN)) +
	geom_boxplot(notch=TRUE) + 
	theme_bw(base_size = 12) + 
	ggtitle("relationship between mRIN (=TIN) and the experimentally determined RIN")

#export the file as PDF landscape
#move all files (tin.xls & plot) back to CAMP server
```

## Assess Nucleotide Content

Converting mRNA to cDNA uses **random primers** which causes certain patterns to be over represented at the beginning (5') end of reads.
We expect A = C = G = T = 25% (if truly random)
![First 10 bases are not equal 25% as expected](http://rseqc.sourceforge.net/_images/NVC_plot.png)
To address this you can trim the first 10 bases of all reads. Then re-perform alignment > assess & compare alignment rates . 

```bash
ml RSeQC
#set BAM input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC

#run read_NVC command on each BAM file using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="read_NVC.py -i $SAMPLE -o ${OUT}/${SRRID}_nucleotide_content"
done
```

## Assess Quality
Use the Phred Quality score to assess quality of base calling. The higher Phred score, the less chance of error.
Quality of bases drops towards the end of the read (3' end) due to synthesis techniques.
Accepted quality is Phred > 30
![Quality distribution](http://rseqc.sourceforge.net/_images/36mer.qual.plot.png)

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/quality

#run each BAM file into read_quality using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="read_quality.py -i $SAMPLE -o ${OUT}_${SRRID}"
done
```
Can trim bases with phred <30.

## Assess PCR duplication

Duplicate reads are reads with the same start/end positions & exact same sequence. You dont expect reads to pile up with the same start & end - this suggests it is from PCR. Whilst you want to remove PCR duplicates you dont want to remove unique reads.

Is it transcription or PCR? There are no fixed start points for DNA replication so with DNA-seq it is well accepted to collapse duplicates. Since there are fixed start points for transcription it is common to see different (non duplicate) reads with the same start sequence.

![PCR duplication](http://rseqc.sourceforge.net/_images/duplicate.png)
A good sample will have the inflexion point low down (i.e. less reads with duplications). The higher the inflexion, the worse the sample.

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/duplication

#run each BAM file into read_duplication using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="read_duplication.py -i $SAMPLE -o ${OUT}_${SRRID}"
done
```

## Assess Sequencing Depth

In DNA-seq easily determine sequencing depth by assessing average coverage over the sequenced region.
In RNA-seq its challenging due to variability in transcript abundance. Use splice junctions detection rate to identify desired sequencing depth. Use a saturation analysis to assess this - check for splcing junctions by sampling 5% of reads, then 10%, then 15% etc until you dont see any new splicing junctions identified (i.e. saturated).
Crucial to ensure there is adequate depth to perform alternative splicing analysis.
![enter image description here](http://rseqc.sourceforge.net/_images/junction_saturation.png)
```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/sequencing_depth

#run each BAM file into junction_saturation using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="junction_saturation.py -i $SAMPLE -r $BED -o ${OUT}_${SRRID}"
done
```

## Assess base distribution

What is the composition of aligned bases? Are they coding (exon) or noncoding. mRNA reads should mostly overlap with exons. The read_distribution.py script calculates how mapped reads are distributed over genome features:
- CDS exon
- 5' UTR exon
- 3' UTR exon
- intron
- intergenic regions

The distribution will depend on the library preparation protocol selected.

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed

#run each BAM file into read_distribution using a For Loop
for SAMPLE in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="read_distribution.py  -i $SAMPLE -r $BED"
done

#print the output to the terminal
more slurm-**********
#Output is printed into the terminal
Total Reads                   10952306
Total Tags                    12476273
Total Assigned Tags           11820580
=====================================================================
Group               Total_bases         Tag_count           Tags/Kb
CDS_Exons           38890108            6051962             155.62
5`UTR_Exons         40728595            427657              10.50
3'UTR_Exons         69257701            2703660             39.04
Introns             1556956010          2520202             1.62
TSS_up_1kb          32291585            8727                0.27
TSS_up_5kb          143148651           29969               0.21
TSS_up_10kb         254502721           45395               0.18
TES_down_1kb        34286667            20436               0.60
TES_down_5kb        147689423           52724               0.36
TES_down_10kb       258048637           71704               0.28
=====================================================================
```

## Assess Insert Size
This is important for paired end reads only. The inner mate (or insert) = the segment of the transcript between the paired end reads R1 & R2.
http://thegenomefactory.blogspot.com/2013/08/paired-end-read-confusion-library.html
```html
fragment                  ========================================
fragment + adaptors    ~~~========================================~~~
SE read                   --------->
PE reads                R1--------->                    <---------R2
unknown gap                         ....................

PE reads      R1--------->                    <---------R2
fragment     ~~~========================================~~~
insert          ========================================
inner mate                ....................
```
Insert size distribution plot: ideally you want a normal distribution of inserts. A bimodal distribution suggests you are selecting another fragment size.
![enter image description here](http://rseqc.sourceforge.net/_images/inner_distance.png)
Inner distance = insert size.
The minus -50 size size happens because there is overlap (the reads R1 & R2 are running into each other). Ideally you want some insert size in between:
```html
fragment          ~~~========================================~~~
insert               ========================================
R1                   ------------------------->                    
R2                                   <-----------------------
overlap                              ::::::::::
stitched SE read     --------------------------------------->
```

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/insert_size

#run each BAM file into inner_distance using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="inner_distance.py -i $SAMPLE -o ${OUT}_${SRRID} -r $BED"
done
```

###  Basic Alignment stats bam_stat.py

ml SAMtools
ml RSeQC

Convert RSeQC Bam_stat.py & samtools flagstat reports into txt file:
```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
#set output directory as alignment QC folder where MultiQC will be run
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC

# run RSeQC bam_stat.py & samtools flagstat commands on each BAM file using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=40GB --wrap="bam_stat.py -i $SAMPLE > ${OUT}/${SRRID}_bam_stat.txt"
	sbatch -N 1 -c 8 --mem=40GB --wrap="samtools flagstat $SAMPLE > ${OUT}/${SRRID}_flagstat.txt"
done
```
These files will then be incorporated into MultiQC report.

# QoRTs (Quality of RNA Seq Toolset)
ml QoRTs

`QoRTs` is a `jar` software file that is an alternative to `RSeQC` that provides a comprehensive & multifunctional toolset assess quality control & data processing of high throughput RNA-seq. It also performs read abundance counts (see Read Quantification chapter). 

QoRTs package is composed of 2 parts: java jar-file (for data processing) & R package (for generating tables, figures, plots)

You will get a return from terminal:
`To execute the QoRTs JAR run: java -jar $EBROOTQORTS/QoRTs.jar` 
Copy and paste the `java -jar $EBROOTQORTS/QoRTs.jar` and place before the QoRTs command. 
[Helpfile](https://hartleys.github.io/QoRTs/jarHtml/index.html) here.
http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf
Also see the [walkthrough example](http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf).

## QoRTs QC command 

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set GTF reference annotation
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
#set output directory
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC

#run QoRTs command for each BAM file using a For Loop
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="java -jar $EBROOTQORTS/QoRTs.jar QC --generatePlots --singleEnded $SAMPLE $GTF ${OUT}/${SRRID}_QoRTs"
done
```
- assumes data is paired unless include `--singleEnded`
- Can run individual functions by specifying their names eg `--runFunctions writeGeneBody` runs only the genebody coverage function. To add further individual functions use a comma without space (comma-delimited list).
- Can exclude individual functions eg `--skipFunctions JunctionCalcs` will run all functions except JunctionCalcs

This is an [example](http://chagall.med.cornell.edu/RNASEQcourse/QC.multiPlot.pdf) of QoRTs output.

### Write QoRTs Decoder Text File

Write the Decoder file: decoder.by.UID.txt
• unique.ID: A unique identifier for the row. THIS IS THE ONLY MANDATORY FIELD. eg SRR5483788
• lane.ID: The ID of the lane or batch. By default this will be set to ”UNKNOWN”. 
• group.ID: The ID of the ”group”. For example: ”VCP” or ”CTRL”. By default this will be set to ”UNKNOWN”. 
• sample.ID: The ID of the biological sample from which the data originated. Each sample can have multiple rows, representing technical replicates (in which the same sample is sequenced on multiple lanes or runs). By default QoRTs will assume that every row comes from a separate sample, and will thus set the sample.ID to equal the unique.ID. 
• qc.data.dir : The directory in which the java utility is to save all the QC data. If this column does not exist, by default it will be set to the unique.ID. 
• input.read.pair.count: The number of reads in the original fastq file, prior to alignment. Find this in the STAR output Log.final.out file.
• multi.mapped.read.pair.count: The number of reads that were multi-mapped by the aligner.

Example:
```
unique.ID group.ID	qc.data.dir  input.read.pair.count
SRR5483788  VCP QoRTs_SRR5483788  12475620
SRR5483789  VCP QoRTs_SRR5483789  16605828
SRR5483790  VCP QoRTs_SRR5483790  13406845
SRR5483794  CTRL  QoRTs_SRR5483794  8038427
SRR5483795  CTRL  QoRTs_SRR5483795  7822813
SRR5483796  CTRL  QoRTs_SRR5483796  9226252
```

### Run QoRTs in Rstudio
```r
library(QoRTs)

#Read in the QC data: 
res <- read.qc.results.data("/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/", 
			decoder.files = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.byUID.txt", 
			calc.DESeq2 = TRUE, calc.edgeR = TRUE); 

# Once you have read in the QC data, you can build all sorts of plots.
# EXAMPLE 1: The makeMultiPlot.all can be used to automatically generate a full battery of multi-plot figures: 
makeMultiPlot.all(res, 
				outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPlots/", 
				plot.device.name = "png"); 

#EXAMPLE 2: Some users may find the large png files difficult to read. QoRTs offers multi-page pdf reports as an alternative, simply by using the plot.device.name parameter: 
makeMultiPlot.all(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPDFs/",
				 plot.device.name = "pdf"); 

#EXAMPLE 3: To print all the basic plots as seperate pngs, use the command: 
makeMultiPlot.basic(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/basicPlots/", 
					separatePlots = TRUE);
					
# Extract size factors. QoRTs generates these to normalise all samples to a comparable scale allowing downstream comparison with DESeq2 or edgeR
get.size.factors(res, outfile = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/sizeFactors.GEO.txt");
```
This also produces count files: QC.geneCounts.txt.gz. It also produces files formatted for DESeq

# MultiQC
ml MultiQC

 Generate a comprehensive & interactive report of post-alignment QC using MultiQC from different Quality Control tools eg RSeqQC, QoRTs. [MultiQC](http://multiqc.info) aggregates results from bioinformatic analyses across samples into a single report
 
 MultiQC searches a folder for analyses & compiles a HTLM report that summarises the output from multiple bioinformatic tools

General post-alignment QC:
- STAR log files
- samtools flagstat
- RSeQCs bam_stat.py

RNA specific QC:
- read distribution (RSeQC or QoRTs)
- gene body coverage (RSeQC or QoRTs)
- splice junction info obtained with QoRTs

Go to the `alignment_QC` folder with the aligned QC files in and run: `multiqc .`

Interpret the [HTML report](https://www.youtube.com/watch?v=qPbIlO_KWN0).
 Compare the post alignment MultiQC HTML reports (the raw unprocessed aligned read report & the trimmed, filtered & depleted aligned read report)

Alternatively to visualise the output of multiple RSeQC reads download the relevant txt files and follow this [R script](https://github.com/friedue/course_RNA-seq2015/blob/master/02_Alignment_QC_visualizeReadDistributionsAsBarChart.R).

# Alignment Assessments
 
Check that alignment rate of RNA-seq reads is > 70%:
 
 Check the aligner's output: `cat FILENAME_Log.final.out`
 - most important number = **uniquely mapped reads**
 - if using >2 BAM files then visualise alignment rate for each using MultiQC

Mount files onto laptop: Right click on Finder --> Connect to server --> Connect to Luscombe Lab
- `infiles <- list.files(path="/Volumes/lab-luscomben/working/oliver/projects/rna_seq_worksheet/alignment_STAR", pattern = "WT_1_Log.final.out", full.names = TRUE)`
- `align.results <- lapply(infiles, function(x) read.table(x, sep="|", strip.white = TRUE, stringsAsFactors = FALSE, skip = 3, fill = TRUE, header = FALSE))`
- `typeof(align.results)`
- `head(align.results[[1]])`
-  `align.results <- lapply(align.results, function(x) transform(x,V2 = as.numeric(gsub("%", "", x$V2) )))`


### Visualise STAR alignment information in R studio using ggplot2
Full explanation available [here](https://github.com/friedue/course_RNA-seq2015/blob/master/01_Alignment_visualizeSTARresults.pdf).

**ggplot2** https://ggplot2.tidyverse.org/
- create graphics : provide data --> map variables to aesthetics --> choose graphical primitives.

R studio script `rna_seq_worksheet_rstudio_script.R` is saved in /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet

Compare the results of STAR alignment across samples:

1. Mount STAR files onto laptop: Right click on Finder --> Connect to server --> Connect to Luscombe Lab

2. Assign variables (run the script):
- plotting correlation
- plotting alignement results
- extract legend
- extract histo info

3. Read in a STAR log file
4. Combine the STAR log files

5. Make data suitable for ggplot2 plotting
- remove % symbols to keep just numbers.
- change appearance of each data frames name
- concatenate (link together in a chain) data frames of align.results
- remove lines without values
- add columns with info on sample & replicate ID using info from row names in the align.results file
- check results had 4 columns

6. Combine various ggplot2 figures into 1
- define the variables that we want to include
- generate bar chart


# Check Mutation Status of Samples (Variant Calling)
[https://www.biostarhandbook.com/introduction-to-variant-calling.html](https://www.biostarhandbook.com/introduction-to-variant-calling.html)
http://oliverelliott.org/article/bioinformatics/wik_bioinform/

Variant calling= identify differences between the observed sequencing reads and a reference genome.

Using the alignment .bam files determine variants at specified genomic locations of interest.
Definitions:
Ploidy = no. of complete chromosome sets in a cell = number of alleles (variant forms).
Genotype = assign one or more reads to a known groups of variants.

**Tools for calling variants:**
- bcftools
- FreeBayes
- GATK (genome analysis tool kit) [https://software.broadinstitute.org/gatk/](https://software.broadinstitute.org/gatk/)
- VarScan2

## VCF: Variance Call Format
https://www.biostarhandbook.com/vcf.html#vcf
VCF Poster:  [http://vcftools.sourceforge.net/VCF-poster.pdf](http://vcftools.sourceforge.net/VCF-poster.pdf)
VCF short summary:  [http://www.htslib.org/doc/vcf.html](http://www.htslib.org/doc/vcf.html)
VCF Specification:  [http://samtools.github.io/hts-specs/](http://samtools.github.io/hts-specs/)
A detailed description of the VCF:  [What is a VCF file](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)  on the GATK forums.

Standard representation for variants is VCF format. All variant callers produce VCFs from BAM alignment files. Samples are represented as a column. A VCF has header & record sections. It is a plain text file but in a specific order.
Can view VCF files in IGV: [http://software.broadinstitute.org/software/igv/viewing_vcf_files](http://software.broadinstitute.org/software/igv/viewing_vcf_files)

![enter image description here](https://lh3.googleusercontent.com/BAGxbWH7dXvbsWulTWG31B7o_SEwt6IAb3UzWRTKKaNTUQV-JZ6mcqPSLMS-8XkM9ICGICq-1s-H5A)

**VCF Header**
Located at the beginning of the VCF file and consists of lines starting with `##`:
```bash
##fileformat=VCFv4.1
##contig=<ID=AF086833,length=18959>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```
The header defines all the acronyms used in the VCF e.g. AF = Allele Frequency.

**VCF Record**
The first nine columns of a VCF file are:

1.  `CHROM`: The chromosome (contig) on which the variant occurs
2.  `POS`: The genomic coordinates on which the variant occurs. For deletions, the position given here are on of the bases preceding the event.
3.  `ID`: An identifier for the variant (if it exists). Typically a dbSNP database if that is known.
4.  `REF`: The reference allele on the forward strand at position `POS`.
5.  `ALT`: The alternate allele(s) on the forward strand. More than one may be present.
6.  `QUAL`: A probability that the  `REF/ALT`  variant exists at this site. It is in Phred scale, just as the FASTQ quality and the MAPQ field in the SAM file are.
7.  `FILTER`: The name of filters that the variant fails to pass, or the value  `PASS`  if the variant passed all filters. If the  `FILTER`  value is  `.`, then no filtering has been applied to the record.
8.  `INFO`: Contains the site-specific annotations represented as  `ID=VALUE`  format. **Represents all samples included in VCF - gives 1 value per VCF.**

*INFO/AD* .. Total allelic depth (Number=R,Type=Integer)
*INFO/ADF* .. Total allelic depths on the forward strand (Number=R,Type=Integer)
*INFO/ADR* .. Total allelic depths on the reverse strand (Number=R,Type=Integer)

9.  `FORMAT`: **Represents individual samples within VCF - gives 1 value per sample.** Sample-level annotations as colon separated TAGS:

This field specifies the meaning of the numbers in each sample column. The TAGS are colon  `:`  separated and map each field of the  `FORMAT`  to each value in the sample column. Suppose that columns 9,10 and 11 of a VCF file were:

```bash
FORMAT         sample1         sample2
GT:PL         0/1:51,0,48    1/1:34,32,0    
```
The variant observed for  `sample1`  has the values  `GT=0/1`  and  `PL=51,0,48`. This same variant when observed in  `sample2`  has the values  `GT=1/1`  and  `PL=34,32,0`
What if you wanted to know what do  `GT`  and  `PL`  mean? Go to the header for their definition:

```bash
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
```

*FORMAT/AD* .. Allelic depth (Number=R,Type=Integer)
*FORMAT/ADF* .. Allelic depths on the forward strand (Number=R,Type=Integer)
*FORMAT/ADR* .. Allelic depths on the reverse strand (Number=R,Type=Integer)
*FORMAT/DP* .. Number of high-quality bases (Number=1,Type=Integer)
*FORMAT/SP* .. Phred-scaled strand bias P-value (Number=1,Type=Integer)

*FORMAT/DV* .. Deprecated in favor of FORMAT/AD;
        Number of high-quality non-reference bases, (Number=1,Type=Integer)
*FORMAT/DP4* .. Deprecated in favor of FORMAT/ADF and FORMAT/ADR;
        Number of high-quality ref-forward, ref-reverse,
        alt-forward and alt-reverse bases (Number=4,Type=Integer)
*FORMAT/DPR* .. Deprecated in favor of FORMAT/AD;
        Number of high-quality bases for each observed allele (Number=R,Type=Integer)
*INFO/DPR* .. Deprecated in favor of INFO/AD;
        Number of high-quality bases for each observed allele (Number=R,Type=Integer)
        
Whereas the  `REF`  and  `ALT`  columns show the change, we need to know how many of the copies of the DNA carry the variant or variants. The variant is encoded in `GT`  that indicates the genotype of this sample at this site. It is constructed out of slash-separated numbers (0/0, 0/1, 1/1) where:
-   `0`  indicates the  `REF`  field,
-   `1`  indicates the first entry in the  `ALT`  field,
-   `2`  indicates the second entry in the  `ALT`  field and so on.

For example for a diploid organism the GT field indicates the two alleles carried by the sample:
-   `0/0`  - the sample is a homozygous reference
-   `0/1`  - the sample is heterozygous, carrying one of each the REF and ALT alleles
-   `1/2`  - would indicate a heterozygous carrying one copy of each of the ALT alleles.
-   `1/1`  - the sample is homozygous for the first alternate
-   `2/2`  - the sample is heterozygous for the second alternate

and so on. For samples of higher ploidy (multiple copies of chromosomes), there would be as many  `/`characters as there are copies of each chromosome. For example, a tetraploid genotype could be  `0/1/0/2`.

 Further columns from 10 onwards represent samples. A VCF file may contain any number of sample columns, thousands even and can be thought of as a single database that represents all variations in all samples.
 
 ## Variant Calling Pipeline

With multiple samples variant call you can either:
1. Call the variants SIMULTANEOUSLY on all BAM files together; or
2. Call variants separately on each BAM file and then MERGE the VCFs

Here we stick to the first approach as merging VCFs is error pone.
```bash
# setting up environments, including paths
ml SAMtools
ml BamTools

FASTA=~/working/oliver/genomes/index/UCSC/hg38.fa
BAM=~/working/oliver/projects/vcp_fractionation/alignment/roi_bam/*.bam
OUT=~/working/oliver/projects/vcp_fractionation/alignment/vcp_mutation_analysis/bcftools

#################
# import region of interest BAM files from sango with rsync
rsync -aP sango-ext:/work/LuscombeU/VCP_nuc_cyto_study/ROI/ .

# remove last two files starting {out} - unclear what these are from
rm {out}ID519_A1_CTRL1_iPS-D0-Cyto_L001.bam {out}ID519_A1_CTRL1_iPS-D0-Cyto_L001.bam.bai

######### samtools mpileup > bcftools view #######

### run as merged samples into 1 VCF
# R155C loc 35,065,364. GLIA, GLIB, CBID. Mutation G > A (RNA mutation is C > T). https://www.ncbi.nlm.nih.gov/clinvar/variation/8469/
samtools mpileup -D -S -d 5000 -r chr9:35065364-35065364 -uf $FASTA *.bam | bcftools view > R155Cstatus_all_samples_merged.vcf

# R191Q genomic location = 35,065,255. CB1E. RNA Mutation C > T (DNA mut is G > A). https://www.ncbi.nlm.nih.gov/clinvar/variation/8473/
samtools mpileup -D -S -d 5000 -r chr9:35065255-35065255 -uf $FASTA *.bam | bcftools view > R191Qstatus_all_samples_merged.vcf

###############
### Extract info from VCF
# check genomic location matches file name
bcftools query -f '%CHROM %POS %REF %ALT\n' R155Cstatus_all_samples_merged.vcf
bcftools query -f '%CHROM %POS %REF %ALT\n' R191Qstatus_all_samples_merged.vcf
# list samples stored in VCF
bcftools query -l R155Cstatus_all_samples_merged.vcf
bcftools query -l R191Qstatus_all_samples_merged.vcf
# extract per sample information
## extract PL flag. GT flag not present - unclear why but not necessary to have as order of PL flag indicates Genotype. Transpose to make human readable.
bcftools query -H -f '%CHROM\t%POS[\t%PL\t]\n' R155Cstatus_all_samples_merged.vcf | rowsToCols -varCol stdin stdout > R155Cstatus_PLflag.vcf
bcftools query -H -f '%CHROM\t%POS[\t%PL\t]\n' R191Qstatus_all_samples_merged.vcf | rowsToCols -varCol stdin stdout > R191Qstatus_PLflag.vcf

#### SAMPLES TO CHECK IN IGV:
ID519_A11_CTRL1_electrically-active-MNs-D35-Cyto_S11_L001.bam
ID519_A11_CTRL1_electrically-active-MNs-D35-Cyto_S11_L002.bam
ID519_E3_GLIA_D3-Cyto_L002.bam
```
## Plot allele frequencies in VCF files
[https://github.com/sndrtj/afplot](https://github.com/sndrtj/afplot)

Use `afplot regions` to plot the single region as a histogram or scatter plot.

Provide on VCF per command > gives one plot output.

VCFs need:
1. AD in FORMAT field
2. indexing with tabix]

R155C_all_samples.vcf
R191Q_all_samples.vcf

####  prepare .txt output from VCF with allele frequencies 



#### use .txt to plot the QC output md file

```




<!--stackedit_data:
eyJoaXN0b3J5IjpbLTU2MTYyOTQ5MiwtMTI3MTgyOTcxNSw0OT
I0NDU5NDIsMjAzMTYxODMwOSwtMTE1MDgyNTA5Nyw2ODE5MTEw
MzMsMTE3MDcyMjMxNCwtNTk4NjcwMzg5LC0xNjExNDYyMTM5LC
0xODYxNzIxOTc3LDExNTkwMTI2MDMsNjc0Mzg4NDM0LC02NDAy
NzE2MzcsODE5MzI5Mzg2LDU2MTkwNDg3NiwtOTIyNjc2NjI2LD
ExNDY4MjczNCwyMDU1MDM5MTY1LDE4NzU3NjkxMTMsMTYwNjgw
NDcwN119
-->