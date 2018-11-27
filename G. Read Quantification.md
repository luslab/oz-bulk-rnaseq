> # Read Quantification
Whilst IGV is used to initial glance at coverage, these methods normalise & objectively quantify gene expression.  To compare the expression rates of individual genes between samples you need to quantify the number of reads per gene.  We first use RNA seq to determine the **abundance** of mRNA (cDNA) fragments, rather than the composition of the fragments. 

2 different ways to quantify mRNA abundances of known genes and transcripts:
1.  Raw Counts - simply the number of reads overlapping with a transcript.
2.  Normalise for gene length & sequencing depth:
	- FPKM (RPKM): Reads/Fragments per kilobase of transcript per millions of read mapped. (FLAWED MEASURE STATISTICALLY). This has been superseded by:
	- TPM: Transcripts per Million: also normalises sequencing depth & gene length but in the reverse order. 

https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# Raw Counts

Instead of normalising simply assign fragments to a defined set of genes/transcripts & determine raw counts. Does not have a step to calculate assembly and different isoforms.

Use a BAM file & GTF file and assign each read as best as possible to a known gene to calculate counts. Then run statistical methods on these counts for differental expression.

Tool = featureCounts, QoRTs (Nobby uses this as it also does QC simultaneously), STAR (also does counts if you give it GTF file), [HTSeq count](http://htseq.readthedocs.io/en/release_0.10.0/index.html).

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

## QoRTs
ml QoRTs

see Chapter 10.2 on Merging Count Data in http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf
see Chapter 7 (page 12) of [walkthrough example](http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf)

QoRTs package is composed of 2 parts: java jar-file (for data processing) & R package (for generating tables, figures, plots)

Write the Decoder file: decoder.by.UID.txt (identical to that used in QC of Aligned Reads chapter)

If there are technical replicates then merge at this point. QoRTs allows count data to be combined across technical replicates. See step 4 (chapter 9, page 15) http://hartleys.github.io/QoRTs/doc/example-walkthrough.pdf 
If there are no technical replicates (as with Nat Comms paper) then skip this step.

```bash
#set QoRTS QC input
QC=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/QoRTs_SRR5483788/,./QoRTs_SRR5483789/,./QoRTs_SRR5483790/,./QoRTs_SRR5483794/,./QoRTs_SRR5483795/,./QoRTs_SRR5483796/
#set output directory
OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/QoRTs_counts

#run QoRTs command for each QC file
java -jar $EBROOTQORTS/QoRTs.jar mergeCounts --mergeFiles $QC --verbose $OUT
```


```r
library(QoRTs)
suppressPackageStartupMessages(library(DESeq2))

decoder.bySample <- read.table( 
					"/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.bySample.txt", 
					header=T,stringsAsFactors=F); 

directory <- "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/"; 
sampleFiles <- paste0(
					 decoder.bySample$sample.ID,
					  "/QC.geneCounts.formatted.for.DESeq.txt.gz" ); 

sampleCondition <- decoder.bySample$group.ID; 
sampleName <- decoder.bySample$sample.ID; 
sampleTable <- data.frame(sampleName = sampleName, 
							fileName = sampleFiles, 
							condition = sampleCondition); 

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
									directory = directory, 
									design = ~ condition); 
dds

dds <- DESeq(dds);
res <- results(dds);
res;

write.table(res, file = "outputData/analyses/DESeq/DESeq.results.txt");

```






## HTSeq-Count
ml HTSeq
ml Pysam

HTSeq-Counts is slow as you cant multithread.

-   [http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

```bash
mkdir -p htseq
#set GTF annoation
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
#set BAM input file
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
#set output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/htseq/htseq_counts

for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 8 --mem=24GB --wrap="htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $SAMPLE $GTF > ${OUT}_${SRRID}.tsv"
done

# merge results files into a single matrix for use in edgeR
join SRR5483788.tsv SRR5483789.tsv | join - SRR5483790.tsv | join - SRR5483794.tsv | join - SRR5483795.tsv | join - SRR5483796.tsv > htseq_counts_table_temp.tsv
echo "GeneID SRR5483788 SRR5483789 SRR5483790 SRR5483794 SRR5483795 SRR5483796" > header.txt
cat header.txt htseq_counts_table_temp.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > htseq_counts_table.tsv
rm -f htseq_counts_table_temp.tsv header.txt
head htseq_counts_table.tsv
```

This output is then analysed for differential expression using edgeR (see next chapter)

# Normalise for Gene Length & Sequencing Depth

There are also counting tools that normalise for sequencing depth at this point. However DE analysis (DESeq & edgeR) normalise later so this isn't necessary at this point. 

[TPM has superseded FPKM as it is a fairer statistical measure](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/). 
RPKM Reads per kilobase of transcript per million mapped reads & FPKM Fragments per kilbase of transcript per million mapped reads. Fragment refers to read pairs from paired-end reads (counting fragments and not individual reads )
FPKM normalises for gene size & library depth using: **FPKM = (10^9^ * C ) / (N * L)**
C = number of mappable reads (fragments) for gene/transcript/exon etc
N = total number of mappable reads in the library
L = number of base pairs in the gene/transcript/exon etc (i.e. the size of gene length) 

Tools = [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/papers/) is now outdated and superseded by [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

## Process:
1. Overlapping bundles of fragment are **assembled**
2. Fragments are connected in an **overlap graph**. This graph models variability in fragment count for each gene across replicates (if the fragments are compatible then these are assumed to come from the same genetic locus - it calculates mutually incompatible fragments where splice sites are mid exon or intron)
3. Transcript isoforms are inferred from the **minimum paths** required to cover the graph
4. Abundance of each gene isoform is estimated with a maximum likelihood probabilstic model. The variance estimates of fragment count for each transcript  are statistically tested to report DE genes.

![enter image description here](https://media.nature.com/lw926/nature-assets/nbt/journal/v28/n5/images/nbt.1621-F1.jpg)

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
eyJoaXN0b3J5IjpbLTQyMzgzODY0NCwyMDkxNTcyMTA2LDEyNz
kzODUxMjEsMTMzMjA2MTUyOSwxODk3NDQyNTgwLDE5MTk2MDYw
MTUsMTcxOTMyMDM4NCw1ODk0NDU3MDgsMTU0NjQ0MzcyMiwtNj
MwMTE3MTY4LC01Mzg2MjU4MjUsLTgzMDU4MTA4MywtNzM0NDE1
NDg5LDM2Nzk2MjY4LDQyMzQwMzcwNCwtMzAzMDkxNTgxLC0zOT
Y3NzY4MjYsMTU5MzMzMDgxNiwyMDI3ODM0OTgzLC0xODk4NDg1
MjU4XX0=
-->