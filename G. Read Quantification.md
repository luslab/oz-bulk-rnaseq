> # Read Quantification
We first use RNA seq to determine the abundance of mRNA (cDNA) fragments, rather than the composition of the fragments. 

Different ways to quantify mRNA abundances of known genes and transcripts:
1.  FPKM (RPKM): Reads/Fragments per kilobase of transcript per millions of read mapped.
2.  Raw Counts: The number of reads overlapping with a transcript.

Where as IGV is used as an initial glance at coverage, these methods normalise & objectively quantify gene expression.  - To compare the expression rates of individual genes between samples you need to **quantify the number of reads per gene.**

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

## Cuffmerge

Allows merger of several Cufflinks assemblies together which is required as even with replicates cufflinks will not necessarily assemble the same numbers & structures of transcripts.
It filters out a number of transfrags that are likely artefacts.
Optional to provide a reference GTF to merge novel isoforms & known isoforms to maximse overall assembly quality 

## Cuffdiff

## cummeRbund

Automatically generates many of the commonly used data visualisations including:
- distribution plots
- correlation plots
- MA plots
- Volcano plots
- Clustering, PCA & MDS plots (global relationships between conditions)
- Heatmaps
- gene/transcript level plots showing transcript structures & expression levels 

Input the output from Cufflinks, Cuffmerge & Cuffdiff

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

# Raw Counts

Instead of calculating FPKM simply assign fragments to a defined set of genes/transcripts & determine raw counts.
Does not have a step to calculate assembly and different isoforms.

Use a BAM file & GTF file and assign each read as best as possible to a known gene to calculate counts. Then run statistical methods on these counts for differental expression.

Tool = [HTSeq count](http://htseq.readthedocs.io/en/release_0.10.0/index.html) , DESeq, edgeR

- `featureCounts` counts reads if any overlap is found with a gene. Can exclude multi-overlap reads or include then for each gene that is overlapped. This is a package of Subread so need to `ml Subread` - Biostars advise this.
- `QoRTs` also does counting - Nobby uses this.

## featureCounts Workflow
ml Subread

1. Count reads (estimate abundance) per sample:
```bash

# Create output folder
mkdir -p featureCounts

#set gene coordinates
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf
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

## Spike-in control

- Spike-ins provide a known concentrations of transcripts that we can compare to the experimental samples
- A common spike in product is [ERCC ExFold RNA spike-in control mix](http://data.biostarhandbook.com/rnaseq/ERCC/ERCC-information.pdf) which is added to the experimental samples
- Fold change in transcript expression between 2 samples tells you about the difference between the 2; not about whether they are highly or lowly expressed.
- At lower transcript expression levels accuracy in determining fold change deteriorates. 
<!--stackedit_data:
eyJoaXN0b3J5IjpbMzIwNTI4NDIxLC0yNzk5MjEzODUsMTQzND
U5MDgwMSwtMjA0NTQ0MDY0NSw3MjQ4ODk1MjcsLTE4ODI2MTcw
NjksMTkzMDY3NDE1NiwxODc3NDkzNDQ5LDE2OTMwNDM0MjYsMT
g5NjkwNDQ2NCwtMjAwMjk4NDQ0NSwtMTkyNjkwNjM5MiwxMDA5
MzAwMTQ5LDExNDAzNzA3OTQsLTIwNzAzNjA2MDcsLTE3OTU0MT
UzODIsNjMzOTMwNjA1LC02MTk1NzU4OCw2MDQ2MDA0NTUsMTU3
NDE4ODE1NF19
-->