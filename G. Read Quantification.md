> # Read Quantification

Whilst IGV is used to initial glance at coverage, these methods **normalise** & objectively quantify gene expression.  To compare the expression rates of individual genes between samples you need to quantify the number of reads per gene.  We first use RNA seq to determine the **abundance** of mRNA (cDNA) fragments, rather than the composition of the fragments. 

Counting reads requires us to know the number of reads originating from a transcript, however counting reads is frequently ambiguous as to their transcript of origin because of multi-isoform genes. Precise measurement of transcript counts is not currently possible. Count based methods (DESeq2 / edgeR) study at gene level - they first requre counting reads across genome features using HTSeq or featureCounts which limits informati. 

## Read Count Process
1. Overlapping bundles of reads are **assembled**
2. Reads are connected in an **overlap graph** that models variability in read count for each gene across replicates (if the reads are compatible then these are assumed to come from the same genetic locus - calculates mutually incompatible reads where splice sites are mid exon or intron).
3. Transcript isoforms are inferred from the **minimum paths** required to cover the graph
4. Abundance of each gene isoform is estimated with a maximum likelihood probabilstic model. The variance estimates of reads count for each transcript  are statistically tested to report DE genes.

![enter image description here](https://media.nature.com/lw926/nature-assets/nbt/journal/v28/n5/images/nbt.1621-F1.jpg)
[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/papers/) is now outdated and superseded by [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)


# Raw Counts
Raw referes to non-normalised read counts. Simply assign reads to a defined set of genes/transcripts & determine raw counts. Does not have a step to calculate assembly and different isoforms.

Use a BAM file & GTF file and assign each read as best as possible to a known gene to calculate counts. Then run statistical methods on these counts for differental expression.

**Tools**: 
[HTSeq count](http://htseq.readthedocs.io/en/release_0.10.0/index.html) - Raphaelle uses this. It is slow as you cant multithread.
featureCounts
QoRTs (Nobby uses this as it also does QC simultaneously)
STAR (also does counts if you give it GTF file)

## HTSeq-Count
ml HTSeq
ml Pysam

[https://htseq.readthedocs.io/en/release_0.11.1/](https://htseq.readthedocs.io/en/release_0.11.1/)
```bash
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=24GB --wrap="htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $SAMPLE $GTF > ${OUT}_${SRRID}.tab"
done

### NEW RUNING ALL SAMPLES AT ONCE FROM ALL TIME POINTS
# create one htseq/ output folder in expression/
mkdir -p /home/camp/ziffo/working/oliver/projects/airals/expression/htseq/

#set GTF annoation
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
# set timepoint folders
TIMEPOINT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D*_samples
#set htseq-count output file aligned

## run multiple alignments using in for loop
for SAMPLE in $TIMEPOINT;
do
DAY=`echo $SAMPLE | grep -E -o 'D[0-9]+_samples'`
#set the SAM file to read in (using STAR output)
SAM=$SAMPLE/*.sam
echo "Running timepoint $SAMPLE"
	for REPLICATE in $SAM 
	do
	#set htseq output file
	OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/htseq
	SRRID=`echo $REPLICATE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem 40G --wrap="htseq-count -s reverse $REPLICATE $GTF > $OUT/$DAY_${SRRID}.tab"
	echo "Running sample $REPLICATE"
	done
done
```

### Optional: merge htseq count files into a single matrix for use in edgeR
```bash
join SRR5483788.tsv SRR5483789.tsv | join - SRR5483790.tsv | join - SRR5483794.tsv | join - SRR5483795.tsv | join - SRR5483796.tsv > htseq_counts_table_temp.tsv
echo "GeneID SRR5483788 SRR5483789 SRR5483790 SRR5483794 SRR5483795 SRR5483796" > header.txt
cat header.txt htseq_counts_table_temp.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > htseq_counts_table.tsv
rm -f htseq_counts_table_temp.tsv header.txt
head htseq_counts_table.tsv
```
Raphaelle uses htseq.  She runs htseq on each sample from each time point separately to create one output table per sample. The output of these are then loaded into the R script **SVD_analysis.Rmd** script where we run the Differential Gene Expression & Visualisation steps. This logs read counts, filters out low read counts & normalises read counts. She then runs heirachical clustering, PCA & SVD 

This output can then be analysed for differential expression using edgeR or DESeq2 (see next chapter)

## QoRTs
ml QoRTs

see Chapter 10.2 on Merging Count Data in http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf
see Chapter 7 (page 12) of [walkthrough example](http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf)

QoRTs package is composed of 2 parts: java jar-file (for data processing) & R package (for generating tables, figures, plots).

QoRTs already produced the Counts files during the initial processing QC run and individual files formatted for DESeq. So skip straight to next chapter for DE analysis.

## featureCounts
ml Subread

1. Count reads (estimate abundance) per sample:

```bash
# Create output folder
mkdir -p featureCounts
#set gene coordinates
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
#set BAM input file
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
#set Counts.txt output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/featureCounts/D7_samples/featureCounts/

#run featureCounts on each BAM file separately
for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=24GB --wrap="featureCounts -a $GTF -g gene_name -o ${OUT}_${SRRID} $SAMPLE"
done

#run featureCounts on all BAM files together. By default it uses gene_id in the GTF - override with gene_name
# Using * wildcard list all BAM files into 1 file. The output file contains a column for each sample. 
featureCounts -a $GTF -g gene_name -o $OUT $BAM
```

For each featureCounts command run there are 2 output files:
- featureCounts_results.txt has actual read counts per gene - tab delimited file where the first six columns contain feature specific information and the rest of the columns contain the read counts that overlap with that feature.
- featureCounts_results.txt.sumary gives quick overview of how many reads were assigned to genes. 

2. Find sequences with highest abundance

To find sequences with most hits sort by column 7: `cat counts.txt | sort -rn -k 7 | head`
Output table is in columns as:
```bash
Geneid    Chr   Start   End	  Strand   Length 	 Hits
```



# MultiQC

Run MultiQC in the counts directory to generate a report of expression using MultiQC from different Quality Control tools eg featureCounts, QoRTs. [MultiQC](http://multiqc.info) aggregates results from bioinformatic analyses across samples into a single report
 
 MultiQC searches a folder for analyses & compiles a HTLM report that summarises the output from multiple bioinformatic tools

Go to the `counts` folder with the aligned QC files in and run: `multiqc .`

Interpret the [HTML report](https://www.youtube.com/watch?v=qPbIlO_KWN0).


## StringTie
ml StringTie

StringTie [manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
StringTie is principally designed for Alignment by HISAT2. In STAR when mapping you need to [specify --outSAMstrandField parameter](https://www.biostars.org/p/172393/)

```bash
mkdir -p stringtie
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
GTF_OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/stringtie/
TSV_OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/stringtie/
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam

for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=24GB --wrap="stringtie -p 8 -G $GTF -e -B -o $GTF_OUT_transcripts_$SRRID -A $TSV_OUT_gene_abundances_$SRRID $BAM"
done

#   '-p 8' tells Stringtie to use eight CPUs
#   '-G <known transcripts file>' reference annotation to use for guiding the assembly process (GTF/GFF3)
#   '-e' only estimate the abundance of given reference transcripts (requires -G)
#   '-B' enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
#   '-o' output path/file name for the assembled transcripts GTF (default: stdout)
#   '-A' output path/file name for gene abundance estimates
```
View transcript records & expression values (FPKM):
 ```bash 
 awk '{if ($3=="transcript") print}' $GTF_OUT_${SAMPLE} | cut -f 1,4,9 | less
```
Press 'q' to exit the 'less' display

Create a tidy expression matrix file for StringTie results:
```bash
# cd to stringtie output directory
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl

./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='SRR5483788,SRR5483789, SRR5483790,SRR5483794,SRR5483795,SRR5483796' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='SRR5483788,SRR5483789, SRR5483790,SRR5483794,SRR5483795,SRR5483796' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='SRR5483788,SRR5483789, SRR5483790,SRR5483794,SRR5483795,SRR5483796' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv

head transcript_tpm_all_samples.tsv gene_tpm_all_samples.tsv
```

# Spike-in control

- Spike-ins provide a known concentrations of transcripts that we can compare to the experimental samples
- A common spike in product is [ERCC ExFold RNA spike-in control mix](http://data.biostarhandbook.com/rnaseq/ERCC/ERCC-information.pdf) which is added to the experimental samples
- Fold change in transcript expression between 2 samples tells you about the difference between the 2; not about whether they are highly or lowly expressed.
- At lower transcript expression levels accuracy in determining fold change deteriorates. 

## ERCC expression analysis

Using read counts from FPKM or raw counts, plot the linearity of the ERCC spike-in reads counts versus the known concentration of the ERCC spike-in Mix.

```bash
# 1. Download file with concentrations & fold change differences for the ERCC spike-in reagent
wget http://genomedata.org/rnaseq-tutorial/ERCC_Controls_Analysis.txt
cat ERCC_Controls_Analysis.txt

# 2. Use perl script to organise the ERCC expected values & our obersved counts
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/Tutorial_ERCC_expression.pl
chmod +x Tutorial_ERCC_expression.pl
./Tutorial_ERCC_expression.pl
cat $RNA_HOME/expression/htseq_counts/ercc_read_counts.tsv

# 3. Run R script to produce x-y scatter plot comparing expected and observed values
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/Tutorial_ERCC_expression.R
chmod +x Tutorial_ERCC_expression.R
./Tutorial_ERCC_expression.R ercc_read_counts.tsv
```
To view the resulting figure, navigate to the below URL replacing  **YOUR_IP_ADDRESS** with your IP address:
-   http://**YOUR_IP_ADDRESS**/rnaseq/expression/htseq_counts/Tutorial_ERCC_expression.pdf
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTc1MjAwMjE3NSwtMTc3ODM2NTEwMyw5Nj
c3ODU1MTYsMTk3NjI4ODc5LDE3MzM1OTc1OTEsMTE4NjY5NTc2
NiwxODY0MTA4MjgxLDEwNzUxMjI4ODIsMjgwODk2OTUzLC04Mj
Q5ODMzODgsNTI3Mzc2NzE0LC03MzgzMTM5MzMsMTA3NzQ2NzYy
NywtMjQ1NjEwOTg3LDk4NjMyMDY1OSwtNDIzODM4NjQ0LDIwOT
E1NzIxMDYsMTI3OTM4NTEyMSwxMzMyMDYxNTI5LDE4OTc0NDI1
ODBdfQ==
-->