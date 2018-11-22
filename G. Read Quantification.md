> # Read Quantification
We first use RNA seq to determine the abundance of mRNA (cDNA) fragments, rather than the composition of the fragments. 

Different ways to quantify mRNA abundances of known genes and transcripts:
1.  FPKM (RPKM): Reads/Fragments per kilobase of transcript per millions of read mapped.
2.  Raw Counts: The number of reads overlapping with a transcript.

Where as IGV is used as an initial glance at coverage, these methods normalise & objectively quantify gene expression.


## FPKM
Reads per kilobase of transcript per million mapped reads = Fragments per kilbase of transcript per million mapped reads. Fragment refers to read pairs from paired-end reads (counting fragments and not individual reads )

Transcript expression is proportional to the number of cDNA fragments that originate from it. However:
- fragment number is biased towards larger genes
- total fragment number is related to total library sequencing depth (number of lanes, multiplex reads) 

FPKM normalises for gene size & library depth using:

**FPKM = (10^9^ * C ) / (N * L)**

C = number of mappable reads (fragments) for gene/transcript/exon etc
N = total number of mappable reads in the library
L = number of base pairs in the gene/transcript/exon etc (i.e. the size of gene length) 

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

- [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/papers/)
- [htseq-count](http://htseq.readthedocs.io/en/release_0.10.0/index.html) has 3 modes union, intersection strict, and intersection nonempty (image above). 
- `featureCounts` counts reads if any overlap is found with a gene. Can exclude multi-overlap reads or include then for each gene that is overlapped. This is a package of Subread so need to `ml Subread` - Biostars advise this.
- `QoRTs` also does counting - Nobby uses this.


## Cufflinks

Process:
Overlapping bundles of fragment are **assembled**
Fragments are connected in an overlap graph
Then tries to pin these fragments down to a gene locus in the reference genome. 
Transcript isoforms are inferred from the minimum paths required to cover the graph
abundance of each isoform is estimated with a maximum likelihood probabilstic model

![enter image description here](https://media.nature.com/lw926/nature-assets/nbt/journal/v28/n5/images/nbt.1621-F1.jpg)

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
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTExNzkxNTE4MzksNjMzOTMwNjA1LC02MT
k1NzU4OCw2MDQ2MDA0NTUsMTU3NDE4ODE1NCwtNjQ2Njk0NDk0
XX0=
-->