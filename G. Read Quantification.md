> # Read Quantification
We first use RNA seq to determine the abundance of mRNA (cDNA) fragments, rather than the composition of the fragments. 

2 different ways to quantify mRNA abundances of known genes and transcripts:
1.  Raw Counts - simply the number of reads overlapping with a transcript.
2.  Normalise for gene length & sequencing depth:
	- FPKM (RPKM): Reads/Fragments per kilobase of transcript per millions of read mapped. (FLAWED MEASURE STATISTICALLY). This has been superseded by:
	- TPM: Transcripts per Million: also normalises sequencing depth & gene length but in the reverse order. 

https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

Whilst IGV is used to initial glance at coverage, these methods normalise & objectively quantify gene expression.  To compare the expression rates of individual genes between samples you need to quantify the number of reads per gene.

# Raw Counts

Instead of calculating FPKM simply assign fragments to a defined set of genes/transcripts & determine raw counts.
Does not have a step to calculate assembly and different isoforms.

Use a BAM file & GTF file and assign each read as best as possible to a known gene to calculate counts. Then run statistical methods on these counts for differental expression.

Tool = [HTSeq count](http://htseq.readthedocs.io/en/release_0.10.0/index.html) , DESeq, edgeR

- `featureCounts` counts reads if any overlap is found with a gene. Can exclude multi-overlap reads or include then for each gene that is overlapped. This is a package of Subread so need to `ml Subread` - Biostars advise this.
- `QoRTs` also does counting - Nobby uses this.

## HTSeq-Count
ml HTSeq

-   [http://www-huber.embl.de/users/anders/HTSeq/doc/count.html](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

```bash
mkdir -p htseq
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
#set BAM input file
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
#set output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/htseq

for SAMPLE in $BAM
do
	sbatch -N 1 -c 8 --mem=24GB --wrap="htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $SAMPLE $GTF > $OUT_${SAMPLE}.tsv"
done

# merge results files into a single matrix for use in edgeR
join SRR5483788_gene.tsv SRR5483789_gene.tsv | join - SRR5483790_gene.tsv | join - SRR5483794_gene.tsv | join - SRR5483795_gene.tsv | join - SRR5483796_gene.tsv > htseq_read_counts_table_all.tsv
echo "GeneID SRR5483788 SRR5483789 SRR5483790 SRR5483794 SRR5483795 SRR5483796" > header.txt
cat header.txt htseq_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > htseq_read_counts_table_all_final.tsv
rm -f htseq_read_counts_table_all.tsv header.txt
head htseq_read_counts_table_all_final.tsv
```

This output is then analysed for differential expression using edgeR (see next chapter)

# FPKM
Reads per kilobase of transcript per million mapped reads = Fragments per kilbase of transcript per million mapped reads. Fragment refers to read pairs from paired-end reads (counting fragments and not individual reads )

Transcript expression is proportional to the number of cDNA fragments that originate from it. However:
- fragment number is biased towards larger genes
- total fragment number is related to total library sequencing depth (number of lanes, multiplex reads) 

FPKM normalises for gene size & library depth using:

**FPKM = (10^9^ * C ) / (N * L)**

C = number of mappable reads (fragments) for gene/transcript/exon etc
N = total number of mappable reads in the library
L = number of base pairs in the gene/transcript/exon etc (i.e. the size of gene length) 

Tool = [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/papers/), Cuffmerge, Cuffdiff, Cuffquants, [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

## featureCounts Workflow
ml Subread

1. Count reads (estimate abundance) per sample:
```bash

# Create output folder
mkdir -p featureCounts

#set gene coordinates
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
#set BAM input file
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/SRR5*_Aligned.sortedByCoord.out.bam
#set Counts.txt output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/featureCounts/D7_samples/counts.txt

#run featureCounts command - by default it uses gene_id in the GTF file. Override this with gene_name attribute.
featureCounts -a $GTF -g gene_name -o counts.txt $OUT $BAM
```
Using the * wildcard you can list all BAM files into 1 text file.

The output file contains a column for each sample. 

For each BAM file there are 2 output files:
- featureCounts_results.txt has actual read counts per gene - tab delimited file where the first six columns contain feature specific information and the rest of the columns contain the read counts that overlap with that feature.
- featureCounts_results.txt.sumary gives quick overview of how many reads were assigned to genes. 

2. Find sequences with highest abundance

To find sequences with most hits sort by column 7: `cat counts.txt | sort -rn -k 7 | head`
Output table is in columns as:
```bash
Geneid    Chr   Start   End	  Strand   Length 	 Hits
```

## Cufflinks
ml Cufflinks

https://github.com/cole-trapnell-lab/cufflinks

Process:
1. Overlapping bundles of fragment are **assembled**
2. Fragments are connected in an **overlap graph**. This graph models variability in fragment count for each gene across replicates (if the fragments are compatible then these are assumed to come from the same genetic locus - it calculates mutually incompatible fragments where splice sites are mid exon or intron)
3. Transcript isoforms are inferred from the **minimum paths** required to cover the graph
4. Abundance of each gene isoform is estimated with a maximum likelihood probabilstic model. The variance estimates of fragment count for each transcript  are statistically tested to report DE genes.

![enter image description here](https://media.nature.com/lw926/nature-assets/nbt/journal/v28/n5/images/nbt.1621-F1.jpg)
```bash
OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/cufflinks/
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf

for SAMPLE in $BAM
do
	sbatch -N 1 -c 8 --mem=24GB --wrap="cufflinks -o $OUT_${SAMPLE} -G $GTF $BAM"
done
```

### Cuffmerge & Cuffdiff

Allows merger of several Cufflinks assemblies together which is required as even with replicates cufflinks will not necessarily assemble the same numbers & structures of transcripts.
It filters out a number of transfrags that are likely artefacts.
Optional to provide a reference GTF to merge novel isoforms & known isoforms to maximse overall assembly quality 


## StringTie
ml StringTie

StringTie [manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
StringTie is designed for Alignment by TopHat and HISAT2. In STAR you need to specify a specific sort.

```bash
mkdir -p stringtie
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
GTF_OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/stringtie/transcripts.gtf
TSV_OUT=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/stringtie/gene_abundances.tsv
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR54837*_Aligned.sortedByCoord.out.bam

for SAMPLE in $BAM
do
	sbatch -N 1 -c 8 --mem=24GB --wrap="stringtie -p 8 -G $GTF -e -B -o $GTF_OUT_${SAMPLE} -A $TSV_OUT_${SAMPLE} $BAM"
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
eyJoaXN0b3J5IjpbLTgwMjE0NzYxOSwyMDQ4MTkwMDQ1LDIxMT
gyNDQzODIsMTEyNTg1MDg0OCwxMTQ4NzE1OTIsLTUzNjE1MTIy
NywtMTIyOTgxNTM3MiwtMTQwNDM3Mzk5MSwtNjYxMDkzMTAwLC
0yNzk5MjEzODUsMTQzNDU5MDgwMSwtMjA0NTQ0MDY0NSw3MjQ4
ODk1MjcsLTE4ODI2MTcwNjksMTkzMDY3NDE1NiwxODc3NDkzND
Q5LDE2OTMwNDM0MjYsMTg5NjkwNDQ2NCwtMjAwMjk4NDQ0NSwt
MTkyNjkwNjM5Ml19
-->