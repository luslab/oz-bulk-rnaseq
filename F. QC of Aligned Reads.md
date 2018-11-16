> # Quality Control of Aligned Reads
 The main STAR output file for QC analysis is **Aligned.sortedByCoord.out.bam**.

After aligning and before performing downstream analyses check for:
1. Excessive amounts of reads not aligned
2. Obvious biases in the read distributions: 3' & 5' bias, base distribution, nucleotide content
3. Read quality
4. Sequencing depth
5. Similarity between replicate samples

All this can be done with RSeQC package http://rseqc.sourceforge.net/

# Bias Assessment
ml RSeQC
ml R

Typical Biases of RNA-seq
- many reads aligned to introns indicates: 
	- incomplete poly(A) enrichment 
	- abundant presence of immature transcripts
- many reads align outside of annoted gene sequences (intergenic reads) indicates:
	- genomic DNA contamination
	- abundant non-coding transcripts
- over representation of 3' portions of transcripts indicates RNA degradation

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

```bash
ml RSeQC
ml R

#set changable elements
## set INDEXED BAM files (BAM.BAI) to read in. can list multiple separated by ","
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/coverage

#run each BAM file into geneBody coverage using a For Loop
for file in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="geneBody_coverage.py -i $file -r $BED -o $OUT -f pdf"
done
```
This produces 2 figures to visualise for 3' or 5' bias. If you detect 3' bias at this stage you can either resequence (costly) or adjust for this bias in downstream analysis.

### Estimate RIN
Colours represent different RIN values (RIN 0 = degraded; RIN 9 = high quality). The RIN 0 line (degraded RNA) shows more 3' bias.
![enter image description here](https://www.researchgate.net/profile/Benjamin_Sigurgeirsson/publication/260841079/figure/fig5/AS:296675668185106@1447744400111/Gene-body-coverage-on-average-for-each-group-Both-RIN-10-and-RiboMinus-show-even.png)

Determine a measure of mRNA degradation in silico using RSeQCs tin.py script to produce a TIN.
- TIN 0 (worst) - 100 (best). TIN 60 = 60% of transcript has been covered.
- tin.py uses the deviation from an expected uniform read distribution across the gene body as a proxy
```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed

#run each BAM file into tin.py using a For Loop
for file in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="tin.py -i $file -r $BED"
done
```
Output is an xls file and a summary txt file (mean & median values across all genes in sample).
Visualise TIN in boxplots in [Rstudio](https://github.com/friedue/course_RNA-seq2015/blob/master/03_mRIN.R) using ggplot

## Assess Nucleotide Content

Converting mRNA to cDNA uses **random primers** which causes certain patterns to be over represented at the beginning (5') end of reads.
We expect A = C = G = T = 25% (if truly random)
![First 10 bases are not equal 25% as expected](http://rseqc.sourceforge.net/_images/NVC_plot.png)
To address this you can trim the first 10 bases of all reads. Then re-perform alignment > assess & compare alignment rates . 

```bash
ml RSeQC
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/*_Aligned.sortedByCoord.out.bam
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/nucleotide_content

#run read_NVC command on each BAM file using a For Loop
for file in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="read_NVC.py -i $file -o $OUT"
done



#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed

#run each BAM file into tin.py using a For Loop
for file in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="tin.py -i $file -r $BED"
done
```



## Assess Quality
Use the Phred Quality score to assess quality of base calling. The higher Phred score, the less chance of error.
Quality of bases drops towards the end of the read (3' end) due to synthesis techniques.
Accepted quality is Phred > 30
![Quality distribution](http://rseqc.sourceforge.net/_images/36mer.qual.plot.png)

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483788_Aligned.sortedByCoord.out.bam,/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483789_Aligned.sortedByCoord.out.bam,/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483790_Aligned.sortedByCoord.out.bam,/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483794_Aligned.sortedByCoord.out.bam,/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483795_Aligned.sortedByCoord.out.bam,/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/SRR5483796_Aligned.sortedByCoord.out.bam
#set designed output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/quality

# run read_quality command
sbatch -N 1 -c 4 --mem=24GB --wrap="read_quality.py -i $BAM -o $OUT"
```
Can trim bases with phred <30.

 ## Assess PCR duplication

Duplicate reads are reads with the same start/end positions & exact same sequence. You dont expect reads to pile up with the same start & end - this suggests it is from PCR. Whilst you want to remove PCR duplicates you dont want to remove unique reads.

Is it transcription or PCR? There are no fixed start points for DNA replication so with DNA-seq it is well accepted to collapse duplicates. Since there are fixed start points for transcription it is common to see different (non duplicate) reads with the same start sequence.

![PCR duplication](http://rseqc.sourceforge.net/_images/duplicate.png)
A good sample will have the inflexion point low down (i.e. less reads with duplications). The higher the inflexion, the worse the sample.

```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/*_Aligned.sortedByCoord.out.bam
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/duplication

# run read_duplication command
read_duplication.py -i $BAM -o $OUT
```

## Assess Sequencing Depth

In DNA-seq easily determine sequencing depth by assessing average coverage over the sequenced region.
In RNA-seq its challenging due to variability in transcript abundance. Use splice junctions detection rate to identify desired sequencing depth. Use a saturation analysis to assess this - check for splcing junctions by sampling 5% of reads, then 10%, then 15% etc until you dont see any new splicing junctions identified (i.e. saturated).
Crucial to ensure there is adequate depth to perform alternative splicing analysis.
![enter image description here](http://rseqc.sourceforge.net/_images/junction_saturation.png)
```bash
#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.bed
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/sequencing_depth

#run junction saturation command
junction_saturation.py -i $BAM -r $BED -o $OUT
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
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed

#run read distribution script
read_distribution.py  -i $BAM -r $BED

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
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/*_Aligned.sortedByCoord.out.bam
#set the reference annotation genome - RSeQC requires BED format (convert GTF > BED)
BED=/home/camp/ziffo/working/oliver/genomes/annotation/Human.GRCh38.GENCODEv24.bed
#set output path & prefix
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/alignment_QC/insert_size

#run inner distance script
inner_distance.py -i $BAM -o $OUT -r $BED
```









### Quality of RNA Seq Toolset (QoRTs)
- `QoRT` is a `jar` software file that is an alternative to `RSeQC` that provides a comprehensive & multifunctional toolset assess quality control & data processing of high throughput RNA-seq.

`ml QoRTs`
You will get a return from terminal like:
`To execute the QoRTs JAR run: java -jar $EBROOTQORTS/QoRTs.jar` 
Copy and paste the `java -jar $EBROOTQORTS/QoRTs.jar` and place before the QoRTs command. 
[Helpfile](https://hartleys.github.io/QoRTs/jarHtml/index.html) here.


`for VARIABLE in *.fastq; do wc -l $VARIABLE; done`

This is a [loop in bash](https://ryanstutorials.net/bash-scripting-tutorial/bash-loops.php). Loops are really helpful to perform a repeated command on multiple things and is built around:
1. provide the list (for X)
2. do the command repeatedly 
3. done when the desired situation is achieved
The `;` provides the sectioning in the command.

_QC command in QoRTs_
`java -jar $EBROOTQORTS/QoRTs.jar QC --singleEnded --seqReadCt 7014609 --generatePdfReport /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam sacCer3.gtf ./QoRTs_output/`

- assumes data is paired unless include `--singleEnded`
- Can run individual functions by specifying their names eg `--runFunctions writeGeneBody` runs only the genebody coverage function. To add further individual functions use a comma without space (comma-delimited list).
- Can exclude individual functions eg `--skipFunctions JunctionCalcs` will run all functions except JunctionCalcs

This is an [example](http://chagall.med.cornell.edu/RNASEQcourse/QC.multiPlot.pdf) of QoRTs output.


# MultiQC

 - Generate a comprehensive report of post-alignment QC using MultiQC from different Quality Control tools eg RSeqQC, QoRTs.
 - [MultiQC](http://multiqc.info) aggregates results from bioinformatic analyses across samples into a single report
 - MultiQC searches a folder for analyses & compiles a HTLM report that summarises the output from multiple bioinformatic tools

General post-alignment QC:
- STAR log files
- samtools flagstat
- RSeQCs bam_stat.py

RNA specific QC:
- read distribution (RSeQC or QoRTs)
- gene body coverage (RSeQC or QoRTs)
- splice junction info obtained with QoRTs

1. collect all QC results of interest into one folder QC_collection
2. create subfolders for each sample
3. run multiQC
`ml multiqc`
`multiqc /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/QC_collection --dirs --ignore ERR* --filename multiqc_align`

Interpret the [HTML report](https://www.youtube.com/watch?v=qPbIlO_KWN0).

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

### Data Simulation
Evaluate performance of processing. 
[Polyester](http://bioconductor.org/packages/release/bioc/vignettes/polyester/inst/doc/polyester.html) is an RNA Seq measurement simulator that simulates fragmentation, reverse-complementing & sequencing.

## ggplot2 R

`ml RSeQC`

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

### Visualise the output of  `RSeQC` quality controls
Basic alignment stats: `bam_stat.py -i WT_1_Aligned.sortedByCoord.out.bam`

To add results of samtools flagstat & RSeQC to a MultiQC report capture the output as a txt file.
`bam_stat.py -i WT_1_Aligned.sortedByCoord.out.bam > bam_stat_WT_1.txt`
`samtools flagstat WT_1_Aligned.sortedByCoord.out.bam > flagstat_WT_1.txt`

To visualise the output of mulple RSeQC reads download the relevant txt files and follow this [R script](https://github.com/friedue/course_RNA-seq2015/blob/master/02_Alignment_QC_visualizeReadDistributionsAsBarChart.R).






<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0MTYyMjQ5MjcsLTQ2MDk3ODQ1MCwtMT
Y0ODMwNjgzMSw2MTMwOTc3MzYsMTQ3NDIxNjQ0MCwyMTA3ODIy
Mzg0LDEzMDM2NzgzMTUsMjIwMzU3Mzk3LDEwMjA2MjU2NjUsMT
E3MTMxMzc5MiwxNzQ3MDk4OTAsLTE4MDkwOTAxNCwtMTE4NDEw
MjA3OCwtMTQ0NDc3NzY2LC0xNDM4MDEzODI5LC0yMTQwMDAxMj
k1LC0xOTcwNDE4OTkwLDYwMzc3MTIwLDE4NDE0NjIxOTgsLTgz
ODE1NDE1OV19
-->