

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
 











 





# Sequence Alignment
- the purpose of alignment (aka pairwise alignment; mapping) is to arrange 2 sequences so the regions of similarity line up. For each base alignment the options are:
	- `-` a space
	- `|` match
	- `.` a mismatch

> Mapping
> -   A mapping is a region where a read sequence is placed.
> -   A mapping is regarded to be correct if it overlaps the true region.
> 
> Alignment
> -   An alignment is the detailed placement of each base in a read.
> -   An alignment is regarded to be correct if each base is placed correctly.

Tools used to be separated into aligners vs mappers. However these have become combined over time. However a dinstiction still exists. For example, studies examining SNPs and variations in a genome would be primarily alignment-oriented. However, studies focusing on RNA-Seq would be essentially mapping-oriented. The reason for this difference is that in RNA-Seq you would need to know where the measurements come from, but you would be less concerned about their proper alignment.

- **CIGAR** (Concise idiosyncratic gapped alignment report string) string is an alignment format used in SAM (sequence alignment map) files. CIGAR uses letters M, I, D etc to indicate how the read aligned to the reference sequence at that specific locus. In  extended CIGAR the symbol `X` is used for mismatches.
**M**  - Alignment (can be a sequence match or mismatch!)
**I**  - Insertion in the read compared to the reference
 **D**  - Deletion in the read compared to the reference
**N**  - Skipped region from the reference. For mRNA-to-genome alignments, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.
**S**  - Soft clipping (clipped sequences are present in read); S may only have H operations between them and the ends of the string
**H**  - Hard clipping (clipped sequences are NOT present in the alignment record); can only be present as the first and/or last operation
**P**  - Padding (silent deletion from padded reference)
 **=**  - Sequence match (not widely used)
**X**  - Sequence mismatch (not widely used)

The sum of lengths of the  **M**,  **I**,  **S**,  **=**,  **X**  operations must equal the length of the read. Here are some examples:
![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/cigar.png)

- Aim of alignment is to identify transcripts (reads) in a sample by mapping them to the genomic origin using a **reference genome**. We want to map millions of reads accurately and quickly. 
- The limitations of alignment are: 
	- sequencing errors
	- genomic variation 
	- repetitive elements
	- multiple different transcript isoforms from same gene
	- main challenge is the spliced alignment of exon-exon spanning reads
	- mapping ambiguity = many reads overlap with more than one isoform
	-  False positives: lowly expressed isoforms are excluded by alignment algorithms —> bias towards identifying strongly expressed genes

![enter image description here](https://www.researchgate.net/profile/Daehwan_Kim13/publication/275410550/figure/fig1/AS:281862078517250@1444212564204/Two-possible-incorrect-alignments-of-spliced-reads-1-A-read-extending-a-few-bases-into.png)


- **Alignment scoring** is based on the value you associate with a match, mismatch or a gap. Alignment algorithms find the arrangement that produce the maximal alignment score. Example scoring matrix:
	- 5 points for a match.
	- -4 points for a mismatch.
	- -10 points for opening a gap.
	- -0.5 points for extending an open gap.
- Modifying the scoring algorithm can dramatically change the way that the sequence is aligned to the reference. See [this example](https://www.biostarhandbook.com/align/misleading-alignments.html) where reducing the gap penalty from -10 to -9 actually corrects the alignment. This setting means that 2 matches (5 + 5 = 10) overcomes a penality gap open (9). Thus to achieve a higher overall score the aligner will prefer opening a gap anytime it can find two matches later on (net score +1). 
- [EDNAFULL](ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4) is the default scoring choice for all aligners.

## Types of alignment

- **Global alignment** is where every base of both sequences has to align to another matching base, to another mismatching base, or to a gap in the other sequence.
- **Local alignment** algorithms look for the highest scoring subregions (or single region). Local alignments are used when we need to find the region of maximal similarity between two sequences.
- **Semi-global alignment** (global-local, global) is a method mixture between the global and local alignments. It attempts to fully align a shorter sequence against a longer one. They are used when matching sequencing reads produced by sequencing instruments against reference genomes. Majority of data analysis protocols rely on semi-global alignment.
- Multiple sequence alignments uses 3 or more sequences (i.e. not pairwise alignment). 
- **Pseudoalignment** - using Kallisto or Salmon - aligns reads to transcripts (not genomes) - much quicker but wont identify new transcripts.

For RNA-seq we need to either :
- align reads to the reference **transcriptome** index (required transcripts to be known and annotated in the reference)
- align reads to the reference **genome** index to identify novel splice events (i.e. reads that cant be aligned to reference transcriptome)

![enter image description here](https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/resize/user/18/Figure19-700x527.png)


### Short read aligners
In 2005 high throughput short-read sequencers changed the face of sequencing which became much cheaper and accessible. Previously sequencing was laborious and expensive and the focus was on producing accurate long reads (1000bp). Sequencing reads longer improves alignment.  Sequencers now produce millions of short reads (50-300bp). Aligners have thus changed to adapt to short reads - rapidly select the best mapping at the expense of not investigating all potential alternative alignments.
- Short read aligners are incredible! They match > 10,000 sequences per second against 3 billion bases of the human genome.
- There is large variation in results between different aligners. A tool that prioritises finding exact matches will do so at the expense of missing locations and overlaps and vice versa.
- Limitations:
	- finds alignments that are reasonably similar but not exact (algorithm will not search beyond a defined matching threshold)
	- cannot handle very long reads or very short reads (<30bp) (become inefficient)

## Splice-aware Aligners

Multiple alignment programmes are available for RNA-seq, each specialising in detecting different factors eg structural variants; fusion transcripts:
-   [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf): Really fast, produces counts for you too. Straight forward RNA seq for differential gene expression analysis. Efficient. Sensitivie. Identified large number of novel splice sites.
-   HISAT2: similar algorithms as Bowtie2/Tophat2 - but performs at much faster speeds
-   Subjunc: designed specifically for gene expression analysis.
-   BWA: the MEM algorithm of bwa supports spliced alignments.
-   BBMap: Is capable of detecting very long splicing.
- Tuxedo suite refers to the pioneering automating pipelines for RNA-seq = alignment with tophat & bowtie, and Quantification & Differential Expression with cufflinks - all now outdated. 
- The new Tuxedo suite refers to
	- alignment: hisat2
	- quantification: stringtie
	- differential expression: ballgown (R package)

When choosing an aligner, we need to decide what features the aligner should provide:
-   Alignment algorithm: global, local, or semi-global?
-   Is there a need to report non-linear arrangements?
-   How will the aligner handle INDELs (insertions/deletions)?
-   Can the aligner skip (or splice) over large regions?
-   Can the aligner filter alignments to suit our needs?
-   Will the aligner find chimeric alignments?

## Reference Genomes
Reference genome sequence repositories:
**[GenCODE](https://www.gencodegenes.org/human/)** - Generally for human data use 
**[ENSEMBL](http://www.ensembl.org/index.html)** : has the most detailed annotations of genomes.
**[UCSC table browser](https://genome.ucsc.edu/)** - Stores large scale result tiles produced by consortia e.g. ENCODE. Contains multiple sequence alignment datasets across entire genomes.https://genome.ucsc.edu/cgi-bin/hgTables
ENCODE
iGenomes
NCBI
Mouse Genome Project
Berkeley Drosphilia Project
RefSeq
[RNA-Central](http://rnacentral.org/) - non coding RNA sequencing database

* UCSC and Ensembl use different naming conventions (which impacts on analyses) - try to stick to one.
* must use Ensembl FASTA file with Ensemble GTF file. You cannot mix Ensembl with UCSC without modifying first and this isn't advised as they have different scaffolds.

### Reference Genome File Formats
Reference sequences are **FASTA files.** Reference sequences are long strings of ATCGN letters.  File formats store start sites, exon, introns. One line per genomic feature.

Reference annotations are **GFF** or **GTF** format.

**GFF = General Feature Format**
9 fields separated by TAB. No white space - all fields except last in each feature must contain value (missing values = `.`)
* Reference sequence: coordinate system of annotation eg Chr1
* Source: annotation software
* Method: annotation type (eg gene)   	[GFF3 = Type: term from lite Sequence Ontology SOFA or accession number]
* Start Position: 1-based integer, <= stop position
* Stop Position: 0-length features (insertion sites)
* Score: sequence identity
* Strand: + = forward. - = reverse. “.” = no stranded
* Phase: codon phase for annotations linked to proteins. 0; 1; or 2 = indicates frame (bases to be removed from beginning of feature to reach first base of next codon)
* Group: class and ID of an annotation which is the parent of the current one.  	[GFF3 = Attributes: TAG=VALUE pairs - ID, name, alias, parent, target, gap]

2 versions in use - similar but not compatible
* version 2 = Sanger Institute http://gmod.org/wiki/GFF2
* version 3 = Sequence Ontology Project http://gmod.org/wiki/GFF3        
 
**GTF = Gene Transfer Format** 
More strict than GFF. Same as GFF2 but 9th field expanded into attributes (like GFF3). http://mblab.wustl.edu/GTF2.html

**BED Format** is the simplest annotation store
3 compulsory fields: chromosome & start & end.
9 optional fields: name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        	field number can vary from 3 - 12. Must be consistent within a file.
indicates region with 0-based start and 1-based end position (GFF & GTF are 1-based in both start and end) Aligning Reads

## Download Reference Sequence (FASTA) & Annotation (GTF) files
Keep reference genomes in a more general location rather than the local folder, since you will reuse them frequently: `/home/camp/ziffo/working/oliver/genomes/sequences` and then make a new directory `mkdir`

**GENCODE process:**
https://www.gencodegenes.org/releases/current.html
find the latest human reference genome: currently this is Release 28 (GRCh38.p12)
Copy the Link address to download the FASTA file to the GENOME SEQUENCE  PRI (primary) - bottom of 2nd table. 
Then in command line download to the appropriate directory: `wget [paste link address]` then `gunzip filename`
Do the same for the GTF file to the **Comprehensive gene annotation PRI (primary) regions**.

 **ENSEMBL process:**
http://www.ensembl.org/info/data/ftp/index.html
Search for species of interest
Click on Gene sets GTF link & DNA FASTA link
GTF: Right click on Saccharomyces)cerevisiae.R64-1-1.92.gtf.gz → copy link address
FASTA: Right click on DNA top level file.
In command line (in appropriate Folder) `wget [paste link address]`
Unzip file `gunzip file_name`

**UCSC process:**
`ml Kent_tools`
Download a GTF file of yeast transcripts from the UCSC Genome Table Browser https://genome.ucsc.edu/cgi-bin/hgTables
Move the downloaded GTF file to the appropriate folder
**![](https://lh4.googleusercontent.com/piQkvTkiSIYCY9m-gATKN8CTmWGFPVZaP7KItC44zJP_oztaNMxjf9O33hljoHvARnSAqaXP1lz5pUo8_7X49xlHKXtX5hUyU-vAfehxNnXAVQ3mh152qUNwlywheUpx5P2GUa4Y)**
N.B. GTF files downloaded from UCSC table have same entries for gene ID and transcript ID → creates problem with analysing different exon isoforms (same gene ID but different transcript ID)

remove first column & first line: `cut -f 2- file_name.txt | sed ‘1d”`
`genePredToGtf file file_name file_name.gtf`

Compress with `gzip` command or `faToTwoBit filename.fa filename.2bit`. Can then reconvert 2bit format —> FASTA format:  `twobittofa file_name.2bit file_name.fa`


## Alignment Workflow:

 - The aim of alignment is to produce a SAM (Sequence Alignment Map) file which contains all information on the sequencing & its alignment. After creating the SAM file there is no need to look at the FASTQ file again since the SAM contains all the FASTQ information.

![RNA-seq Flowchart - Module 1](https://github.com/griffithlab/rnaseq_tutorial/wiki/Images/RNA-seq_Flowchart2.png)

### 1. Build Index

- All short read aligners first build an index from the reference genome. Then the FASTQ sequencing files are aligned against this index.
- Index building prepares the reference genome to allow the tools to search it efficiently. It creates multiple files that should be stored near to the FASTA reference sequence.
- Can take days for large genomes to build an index. The relevant indexes have already been created and stored on the CAMP cluster in `/home/camp/ziffo/working/luscombelab-UCL`
 
Input Files = Reference genome sequence (FASTA) & Annotation file (GTF)

Create directory to store index in: `mkdir GRCh38.12_STAR_index`
Module load STAR `ml STAR`  

STAR commands have the format: `STAR --option1-name option1-value(s)--option2-name option2-value(s) ...`
To generate the index in STAR, specify the location of:
1. index directory to store the output
2. FASTA reference genome file
3. GTF annotation file
4. Overhand: read length minus 1. Read length distribution are shown in the MultiQC report. This is length of the genomic sequence around the annotated junction to be used for the splice junctions database

```bash
#Set the changable elements
IDX=/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf

#Send cmd to generate index as batch job to cluster:
sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runMode genomeGenerate --genomeDir $IDX --genomeFastaFiles $REF  --sjdbGTFfile $GTF --sjdbOverhang 59 --runThreadN 8"
```

The above code is resuable and applicable to all situations by editing the $changable elements.
The "hard-coding" looks like this (very long and difficult to edit):
```bash
sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runMode genomeGenerate --genomeDir  --genomeFastaFiles /home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf --sjdbOverhang 59 --runThreadN 8"
```

Check it is running using `myq`
If it isn’t seen there then `ls` the folder → look for output file named “slurm…”
Open slurm file: `more slurm…` which will explain outcome of file eg FATAL INPUT PARAMETER ERROR

Alternative approach is to build Index using [snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) which parallelises the indexing process: https://www.biostarhandbook.com/unix/makefiles.html

 - Write the Snakefile in Atom and save it in the relevant directory near the index. Define input, output and shell script.
```bash
input:
	ref="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa", starref="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR",
output: 
	"/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR"+"/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index/chrName.txt"
shell: 
	"STAR --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.ref} --sjdbOverhang 59 --runThreadN 8"
```
 
### 2. Align each FASTQ file to the Index
ml Python/3.5.2-foss-2016b
ml SAMtools

* Sample distributed over X flow cell lanes → X fastq files per sample
* STAR merges the X files if multiple file names are indicated (other align tools don't)
* Separate the file names with a comma (no spaces)
* Create directory to store STAR output `mkdir alignment_STAR`

By allocating all file names to the `$INPUT` term it means all the FASTQ files are read together (see example further down) but it is best to do each separately - can use [snakemake](http://slides.com/johanneskoester/snakemake-tutorial#/) - this makes it easier to see if there is an alignment error in each individual sequencing file:

* List fast.qz files separated by commas and remove white spaces:

`INPUT= ls -m /home/camp/ziffo/working/oliver/projects/airals/fastq_files/SRR5*.fastq| sed 's/ //g' | echo $INPUT | sed 's/ //g'`

`ls -m` list, fill width with a comma separated list of entries all the fastq files.
`sed` remove a space from each file name
`echo` displays a line of text - in this case it lists all the file names  - in this case the pipe to echo is to remove new spaces that are created between multiple lines. no space after FILES - with space after it thinks FILES is command.
`sed` = stream editor - modify each line of a file by replacing specified parts of the line. Makes basic text changes to a file - `s/input/output/g`

```bash
#SET CHANAGABLE ELEMENTS
#set the index
IDX=/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index
#set the fastq sequencing file to read in
READ1=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/rRNA_depleted/SRR5483788_1_no_rRNA.fq
#set the paired fastq sequencing file to read in (for paired end data only)
READ2=
#set name under which to store the BAM file output
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/SRR5483788_

#SEND ALIGNMENT AS SBATCH
sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runThreadN 1 --genomeDir $IDX --readFilesIn $READ1 --outFileNamePrefix $BAM --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic"
```
The above script needs editing for each sequencing file (redundant). 
To create a script that runs all sequencing files (non-redundant) you can use a for loop or snakemake.

**For Loop**
```bash
# Create output folder
mkdir -p bam

# Exit this script on any error.
set -euo pipefail

# set the index
IDX=/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index

for SAMPLE in VCP CTRL;
do
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files.
        READ1=reads/${SAMPLE}_${REPLICATE}_R1.fq
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        # Run the aligner.
        STAR --runThreadN 1 --genomeDir $IDX --readFilesIn $READ1 --outFileNamePrefix $BAM --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic
        # Index each BAM file as they are produced
        samtools index $BAM
    done
done

#Place the above into a script with atom & save in appropriate folder, name it `align.sh` and run it with
bash align.sh
```

**Snakefile**
- write the following rules in Atom and save the file in the appropriate directory
```bash
rule star_map:
input: 
	index="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"
	sample="/home/camp/ziffo/working/oliver/projects/airals/fastq_files/{sample}.fastq", 
output: 
	"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam"
shell:
	"STAR --genomeDir {input.index} --readFilesIn {input.sample} --outFileNamePrefix {output}_ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 1"

#dry run workflow:
`snakemake -np /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam`
#execute workflow:
`snakemake /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam`
```
## STAR output

The [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  has all the explanations on how to write the STAR command and fine tune parameters including:
* multi-mapped reads
* optimise for small genomes
* define min & max intron sizes 
* handle genomes with >5000 scaffolds
* detect chimeric & circular transcripts
if files are compressed add `--readFilesCommand gunzip -c`
STAR ENCODE options: `outFilterMultimapNmax 1` max number of multiple alignments allowed for a read - if exceeded then read is considered unmapped i.e. when 1 is used this only identifies unique reads and removes multimapped reads. This is generally accepted.

#`STAR` will perform the alignment, then extract novel junctions which will be inserted into the genome index which will then be used to re-align all reads
#`runThreadN` can be increased if sufficient computational power is available

**STAR Output Files**
 - Aligned.sortedByCoord.out.bam - the loci of each read & sequence
 - Log.final.out - summary of alignment statistics
 - Log.out - commands, parameters, files used
 - Log.progress.out - elapsed time
 - SJ.out.tab - loci where splice junctions were detected & read number overlapping them
 - Unmapped.out.mate1 - fastq file with unmapped reads

1. check the slurm file: `more slurm-132XXX.out`

2. check the summary of the output: `cat FILENAME_Log.final.out`
![summary of alignment](https://user-images.githubusercontent.com/33317454/37438378-095cb442-282d-11e8-95fd-bfecefae5b75.png)

3. check that the number of reads in the fastq files matches the STAR log output: 
`more SRR5483788_Log.final.out`
`echo $( cat /home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/SRR5483789_1.fastq | wc -l)/4 | bc`

### Index BAM read alignments
ml SAMtools

- STAR auto creates [BAM files](http://software.broadinstitute.org/software/igv/bam).  `-b`will produce a BAM file. `-s` will produce a SAM file.
- Now need to index each BAM file. The indexed BAM file format = BAM.BAI file. Inxed allows quick access to the BAM files without having to load them to memory.

With a For Loop above you can automate this with the above phase as each BAM file is produced:
```bash
for file in ~/working/oliver/projects/airals/alignment_STAR/D7_samples/trimmed_filtered_depleted/SRR5*_Aligned.sortedByCoord.out.bam
do
	sbatch -N 1 -c 8 --mem 40 --wrap="samtools index $file";
done
```
for individual files: `samtools index filename_Aligned.sortedByCoord.out.bam`

## Sequence Alignment Maps (SAM) files

[SAM files](https://www.biostarhandbook.com/sam/sam-flags.html) are generic nucleotide alignment format describing alignment of sequenced reads to a reference. SAM format are TAB-delimited line-orientated (each row represents a single read alignment) text consisting of 2 types of tags:
1. Header: meta-data
2. Alignment: longer with information on alignment.

-   The  [SAM format specification](http://samtools.github.io/hts-specs/SAMv1.pdf)  is the official specification of the SAM format.
-   The  [SAM tag specification](http://samtools.github.io/hts-specs/SAMtags.pdf)  is the official specification of the [SAM tags](https://www.biostarhandbook.com/sam/sam-flags.html).
-  Each SAM file has 11 mandatory columns. Despite this aligners vary in how much information they put into these columns - impossible to transform one aligners SAM file into another aligners SAM file.

![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/bam_structure.png)

### SAM header tags
- Header tags start with `@`then value pairs 2 letter abbreviations e.g. `SQ` = sequence directory listing  `SN` sequence name aligned against (from FASTA), `LN`  sequence length, `PG` program version that was ran.
* 1 line per chromosome
* To retrieve only the SAM header `samtools view -H`
* To retrieve both the header & alignment sections `samtools view -h`
* The default `samtools view` will not show the header section

`samtools view -H filename.bam`

### SAM alignment tags
* Each line corresponds to one sequenced read.
* The SAM tags are defined [here](http://samtools.github.io/hts-specs/SAMtags.pdf)
* [11 mandatory columns](https://www.biostarhandbook.com/sam/sam-flags.html) in specified order:

`<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>`

*[2nd column FLAG field](https://www.biostarhandbook.com/sam/sam-flags.html):*
* stores info on the respective read alignment in one single decimal number
* decimal is the sum of all the answers to Yes/No questions:
```
   Binary       Integer   Name              Meaning 
000000000001       1     PAIRED          Read paired
000000000010       2     PROPER_PAIR     Read mapped in proper pair
000000000100       4     UNMAP           Read unmapped
000000001000       8     MUNMAP          Mate unmapped
000000010000      16     REVERSE         Read reverse strand
000000100000      32     MREVERSE        Mate reverse strand
000001000000      64     READ1           First in pair, file 1
000010000000     128     READ2           Second in pair, file 2
000100000000     256     SECONDARY       Not a primary alignment 
001000000000     512     QCFAIL          Read fails platform/vendor quality checks
010000000000    1024     DUP             Read is PCR or optical duplicate
100000000000    2048     SUPPLEMENTARY   Supplementary alignment
```

* To decode the FLAG integer into plain english [click here](https://broadinstitute.github.io/picard/explain-flags.html). 

_Optional fields_
Following the 11 mandatory fields, the optional fields are presented as key-value pairs in the format of  `<TAG>:<TYPE>:<VALUE>`, where  `TYPE`  is one of:
-   `A`  - Character
-   `i`  - Integer
-   `f`  - Float number
-   `Z`  - String
-   `H`  - Hex string

Reads within the same SAM file may have different numbers of optional fields, depending on the aligner program that generated the SAM file. 
`samtools flagstat` assesses the FLAG field and prints a summary report: `samtools flagstat Aligned.sortedByCoord.bam`

### Read Group tags
[Read group](https://www.biostarhandbook.com/sam/sam-analyze.html) tags `RG` contain sample information in the BAM file. e.g. read group ID, library, sample etc.
Print read group tags `samtools view -H filename.bam | cut -f 12-16 | head -1`
You may need to add a custom read group tag in individual bam files prior to merging with `samtools addreplacerg` using `TAG:FORMAT:VALUE`
This allows pooling results of multiple experiments into a single BAM dataset to simplify downstream logistics into 1 dataset.

[BAM readgroups GATK website](https://gatkforums.broadinstitute.org/gatk/discussion/1317/collected-faqs-about-bam-files)

Most important readgroup tags: `ID`, `SM`, `LB`, `PL`

### BAM versus SAM
- BAM file = Binary Alignment Map - human readable TAB-delimited compressed. BAM files are binary, compressed and sorted representation of SAM information with alignment coordinates allowing fast query on info by location. Used to exchange data as it quicker than SAM.
- Bigger than gzipped SAM files as they are optimised for rapid access (not just size reduction). SAM files are human readable, BAM are compressed. BAM are much smaller. 
- for SAM files you can run other commands on them eg head FILENAME.sam whereas BAM files need to be run through samtools i.e. `samtools view FILENAME.bam | cut -f 2 | head`

You can convert BAM > SAM if needed. 
`samtools sort filename_Aligned.sortedByCoord.out.bam > filename_Aligned.sortedByCoord.out.sam`
or alternatively: 
`samtools view -h FILENAME_Aligned.sortedByCoord.out.bam > FILENAME_Aligned.sortedByCoord.out.sam`
Convert a BAM file into a SAM file (including the header): `samtools view -h FILENAME.bam > FILENAME.sam`
Compress a SAM file into BAM format (-Sb = -S -b)" `samtools view -Sb FILENAME.sam > FILENAME.bam`
To peak into a SAM or BAM file: `samtools view FILENAME.bam | head`

### CRAM format
Similar to BAM (binary compressed) but smaller as some compression is in the reference genome.
Sometimes you need the reference genome information so these arnt always appropriate.
Supported by samtools
Concerted effort to move from BAM to CRAM.

## Manipulating BAM files
There are 4 major toolsets for processing SAM/BAM files:

- [SAMTools](http://www.htslib.org/) - interact with high throughput sequencing data, manipulate alignments in SAM/BAM format, sort, merge, index, align in per-position format. SAMTools help page = `samtools --help` Usage:   `samtools <command> [options]`

5 key SAMTool commands:
1. Indexing
2. Editing
3. File operations (aligning, converting, merging)
4. Statistics
5. Viewing
- [Picard](https://broadinstitute.github.io/picard/) - Java tools for manipulating high throughput sequencing data
- [DeepTools](https://deeptools.readthedocs.io/en/develop/) - visualise, quality control, normalise data from deep-sequencing DNA seq
- [BAMtools](https://github.com/pezmaster31/bamtools/wiki/Tutorial_Toolkit_BamTools-1.0.pdf) - read, write & manipulate BAM genome alingment files . `ml BamTools`

### Filter data from BAM files
- Unlike other aligners, STAR already creates separate bam files of aligned and unmapped.
- Remove poor alignments - eg Phred scale quality <20 & keep alignments which are "properly paired"
- Use FLAGS to filter using samtools: `-f` flag includes matches; `-F` flag includes mismatches. Flag `4` indicates unmapped.  `-h` is used to print the header

Count alignments that were aligned: `samtools view -c -F 4 filename.bam`
Get overview of alignments: 
report of flags `samtools flagstat filename.bam` 
report on how many reads align to each chromosome `samtools idxstats filename.bam` 
report on flags `bamtools stats -in filename.bam`

Filter on mapping quality using `-q` flag:
Count number of reads with mapping quality >1 `samtools view -c -q 1 filename.bam`
Create bam file with mapping quality >= 20`samtools view -h -b -q 20 FILENAME.bam > high_mapq_reads.bam` 

Calculate the depth of coverage: `samtools depth filename.bam | head`  or `bamtools coverage -in filename.bam | head`
Count reads that have multiple alignments `samtools view -c -F 4 -f 256 filename.bam`
Count supplementary (chimeric) alignements `samtools view -c -F 4 -f 2048 filename.bam`
Count primary (non supplementary & non secondary) alignments `samtools view -c -F 4 -F 2304 filename.bam`
Create BAM with uniquely aligned reads (STAR gives uniquely aligned reads a mapping quality of 255 so pull reads with mapping quality = 255 `samtools view -h -q 255 FILENAME.bam > uniquely_aligned_reads.bam` .  

Create BAM file with only reads aligned to reverse strand:
- sort the BAM file (A-Z; 0-9) and count the adjacent lines which are identical using `samtools view FILENAME.bam | cut -f 2 | sort | uniq -c`
- Then, create file with the specific feature e.g. reverse reads (FLAG = 16 in column 2): `samtools view -h -f 16 FILENAME.bam > reverse_reads.bam`

Create SAM file with reads of insert sizes > 1000bp use the CIGAR string (column 6 in SAM file):
- first convert BAM to SAM file ` samtools view -h FILENAME.bam > FILENAME.sam`
- use `grep` to exclude (using `-v`) lines with >3 digits (using `[0-9][0-9][0-9][0-9]`) followed by `N` (N means mismatch i.e. skipped bases, whereas M = match) `egrep -v "[0-9][0-9][0-9][0-9]N" FILENAME.sam > smallinsert_reads.sam`
- Alternatively use `awk` to focus on column 6 (`$6` in CIGAR string) and exclude lines with 3 digits (using `![0-9][0-9][0-9][0-9]`) then printing everything `{print $0}` and creating new file:  `awk '!($6 ~ /[0-9][0-9][0-9][0-9]N/) {print $0}' FILENAME.sam > smallinsert_reads.sam` 

Create SAM file with intron spanning reads:  
- use `grep` to select lines with a number of digits (using `[0-9]+`) then `M` (i.e. matches) then any number of digits again, then `N` (i.e. mismatches) then any number of digits and then M again at the end: `egrep "(^@|[0-9]+M[0-9]+N[0-9]+M)" FILENAME.sam > intron-spanning_reads.sam`
- Alternatively use awk to focus on column 6 (CIGAR string) and select the header `$1 ~ /^@/` and the 6th column with any number of digits followed by M followed by digits then N then digits the M: `awk '$1 ~ /^@/ || $6 ~ /[0-9]+M[0-9]+N[0-9]+M/ {print $0}' FILENAME.sam > intron-spanning_reads.sam`

### Merge BAM files

As you aligned each fastq file separately you have a BAM file for each fastq. At some point you will need to merge all the BAM files for downstream processing.  `samtools merge all_bam_files.bam filename1.bam filename2.bam filename3.bam`
Check the new merged bam file: `samtools view -H all_bam_files.bam`

# Data Visualisation

## Genome Browser
- Check results visually to ensure reads align to expected regions without excess mismatches.
- Genome browsers give a linear track representing the forward stand of the genome. left = 5'. right = 3'
- can visualise line orientated formats (fasta, bed, gff, SAM/BAM)
- genomic features are drawn over the linear track in **glyphs** (pictogram)
	- Horizontal intervals = directions, genes, alignments
	- Values over intervals = coverages, probabilities
	- Attributes at locations = mutations, deletions, junctions

#### `samtools tview`
- simplest genome browser. can visualise any DAM file `samtools tview --reference reference_genome.fa filename.bam`

#### Standalone Genome Browsers
- [Integrative Genomics Viewer IGV](http://software.broadinstitute.org/software/igv/book/export/html/6)  by the Broad Institute
- [Integrated Genome Browser IGB](https://bioviz.org/) by University of North Carolina
- [Artemis](https://www.sanger.ac.uk/science/tools/artemis) by the Wellcome Sanger Institute
- [Tablet](https://ics.hutton.ac.uk/tablet/) by James Hutton Institute
- [SeqMonk](https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/) for high throughput RNA-seq
- [iobio](http://iobio.io/) for real-time genomic analysis

#### Online Genome Browsers
- [Ensembl](http://useast.ensembl.org/index.html)
- [UCSC](https://genome.ucsc.edu/)

## [Integrative Genomics Viewer (IGV)](https://www.biostarhandbook.com/visualize/igv.html)
`ml IGVTools`

Best resources are the [IVG mannual](http://software.broadinstitute.org/software/igv/userguide) and [youtube videos](https://www.youtube.com/results?search_query=integrative+genome+viewer)
1. Run IGV on local computer and mount CAMP. [Set Java 8 as default](https://stackoverflow.com/questions/46513639/how-to-downgrade-java-from-9-to-8-on-a-macos-eclipse-is-not-running-with-java-9) since IGV doesnt work with Java 10
On local terminal `cd ~/bin/IGV_2.4.14/lib` & run IGV via command line on local terminal: `java -Xmx750m -jar igv.jar`
2. Set reference genome to Human (hg38) top left box.
3. Click File load from file > click Desktop > mount CAMP locally > click relevant BAM & BAI files (can load multiple at once).

To visualise on IGV its easier to generate TDF files which are much lighter. This is useful if want to add more data-sets later. To generate TDF files first generate Bedgraph coverage files, then sort and then create the tdf file. Create 3 different coverage files: positive, negative strands and total. As the data is stranded it is better to look at both strands separately. Run the code in file: PE_strandedBedGraph.sh

4. Rename BAM files: right click file name in left hand column.
5. Go to the Genomic Location of interest. For SFPQ type chr1:35,176,378-35,193,158 in top middle box > Go. Find Genomic Locations using google search eg https://www.genecards.org/cgi-bin/carddisp.pl?gene=SFPQ
Mark the region of interest: Regions > Region Navigator > Add
6. Zoom in & out using top right zoomer or +/-. IGV will ask you to zoom in to see alignment data. This is because it has a default of 30 kb memory. You can change this to a higher value to see alignment data from a further out zoom but this will slow down the processing speed when scrolling across the genome. If you have deep coverage files then keep memory low to avoid it slowing down processing.

For each BAM file there are 2 tracks. Top track shows summary distribution of the coverage of the exonic islands separated by spice junctions (introns). Bottom tracks show all the individual sequence read alignments piled up.

7. Compare exon usage in the top tracks to transcript genome (bottom row): make the rows for each of the BAM files smaller and right click on the reference genome at the bottom row (in blue) > expand to show all the previously annotated differential splicing.

Bases that dont match the reference sequences are highlighted by colour. The deeper the shade of grey the more confidence you can have that the sequence was aligned correctly. White means no confidence alignment.

![IGV screen layout](https://lh3.googleusercontent.com/sOSIM9tCveT60ZWs2W9-AliTlVfLPmO4ik9w_ZFBAKh90z5HH9qLXRMQQ0RCajk73UL-ypVJYQbw7w "IGV screen layout")

![IGV RNA-seq specific view](https://lh3.googleusercontent.com/h7PbqBtb3kHxxevIpjvKJUAd451K0UFOoACMogIZzUhVVMz-_AqRnjSYsNpmhYeCbct9ikfaZU8-Yg "IGV RNA-seq specific view")
IGV is used only to validate & confirm analysis results.  Use it to explore large genomic datasets. It is not good for the primary analysis.
![IGV of SFPQ D7 NPC samples](https://lh3.googleusercontent.com/r8Ph08oRuLWUBmnc6gbEyX5Rg3iBEkGhNmmNTHqTr7J01dtwdBGIdAqYJ2BMNlLcVIyYxPbn0QEhTQ "IGV of SFPQ D7 NPC samples")

Top row = chromosome. red bar is location. blue lines mid-section refer to transcripts binding with more = higher peak. bottom section = reference genomes.
Coverage line = quick identification of highly covered regions
Grey boxes = aligned reads. 
Gaps between reads with horizontal grey line = introns
Blue lines within reads = insertions
Red lines within reads = deletions
Blue boxes = reference genome

**Sashimi Plots**
Visualise splice junctions & explore exon usage
![sashimi plot explained](http://miso.readthedocs.io/en/fastmiso/_images/sashimi-plot-example-annotated.png)
Bar graph height =  read coverage
Arcs = splice junctions
Numbers = number of reads that contain the respective splice junction.
IGV does not normalise for read number per sample in sashimi plots so dont overinterepret the read counts.

# Quality Control of Aligned Reads
Analyses now switch from command line to R studio.
 The main STAR output file for downstream analyses are **Aligned.sortedByCoord.out.bam** & **Log.final.out**.
**out.mate1 files** are fastx files that are uncompressed and can be made smaller using gzip.

After aligning and before performing downstream analyses check for:
1. Excessive amounts of reads not aligned
2. Obvious biases in the read distributions
3. Similarity between replicate samples

### Alignment Assessments
 
Check that alignment rate of RNA-seq reads is > 70%:
 
 Check the aligner's output: `cat FILENAME_Log.final.out`
 - most important number = **uniquely mapped reads**
 - if using >2 BAM files then visualise alignment rate for each e.g. using MultiQC in **R studio**: see https://github.com/friedue/course_RNA-seq2015/blob/master/01_Alignment_visualizeSTARresults.pdf
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


## Bias Identification

### Typical Biases of RNA-seq
- many reads aligned to introns indicates: 
	- incomplete poly(A) enrichment 
	- abundant presence of immature transcripts
- many reads align outside of annoted gene sequences (intergenic reads) indicates:
	- genomic DNA contamination
	- abundant non-coding transcripts
- over representation of 3' portions of transcripts indicates RNA degradation

### Read distribution
- mRNA reads should mostly overlap with exons. Test this with `read_distribution.py` script
	- counts number of reads overlapping with various genes & transcript associated genomic regions (introns and exons)
- download BED file from [UCSC genome](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=685446505_FqnRnlREChczp8SYDIJOSvLwBshv&clade=other&org=S.+cerevisiae&db=sacCer3&hgta_group=genes&hgta_track=sgdGene&hgta_table=0&hgta_regionType=genome&position=chrIV%3A765966-775965&hgta_outputType=primaryTable&hgta_outFileName=).

`read_distribution.py -r sacCer3.bed -i /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam`

Output: 
Total Reads                   937851
Total Tags                    947664
Total Assigned Tags           0

```
Group             Total_bases         Tag_count           Tags/Kb
CDS_Exons           8832031             0                   0.00
5'UTR_Exons         0                   0                   0.00
3'UTR_Exons         0                   0                   0.00
Introns             69259               0                   0.00
TSS_up_1kb          2421198             0                   0.00
TSS_up_5kb          3225862             0                   0.00
TSS_up_10kb         3377251             0                   0.00
TES_down_1kb        2073978             0                   0.00
TES_down_5kb        3185496             0                   0.00
TES_down_10kb       3386705             0                   0.00
```

Visualise this output using this [R script](https://github.com/friedue/course_RNA-seq2015/blob/master/02_Alignment_QC_visualizeReadDistributionsAsBarChart.R).

### Gene body coverage
Assess 3' or 5' biases using **RSeQC** `geneBody_coverage.py` script:
- uses an annotation file with transcript models of choice
- it divides each transcript into 100 sections
- then counts reads overlapping with each section
- produces 2 plots showing abundance of reads across transcript bodies

generate index for the BAM file:
`samtools index WT_1_Aligned.sortedByCoord.out.bam`

run the script on the aligned reads & annotation file and set the output file name.
`geneBody_coverage.py -i WT_1_Aligned.sortedByCoord.out.bam -r /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/sacCer3.bed -o geneBodyCoverage_WT_1 -f pdf`

CAMP will not produce a png image file. Only a PDF so use `-f pdf` to set output as PDF file.

Produces 2 figures to visualise for 3' or 5' bias:
X axis = gene body percentile (left is 5' end; right is 3' end)
Y axis = coverage
Lines represent different quality RNA (RIN 0 = degraded; RIN 9 = high quality). The RIN 0 line (degraded RNA) shows more 3' bias.

![enter image description here](https://www.researchgate.net/profile/Benjamin_Sigurgeirsson/publication/260841079/figure/fig5/AS:296675668185106@1447744400111/Gene-body-coverage-on-average-for-each-group-Both-RIN-10-and-RiboMinus-show-even.png)

### mRIN calculation using tin.py
- RNA integrity number (RIN) is rarely reported in public data repositories.
- Instead determine a measure of mRNA degradation in silico using RSeQCs tin.py script to produce a TIN.
- TIN 0 (worst) - 100 (best). TIN 60 = 60% of transcript has been covered.
- tin.py uses the deviation from an expected uniform read distribution across the gene body as a proxy

`tin.py -i WT_1_Aligned.sortedByCoord.out.bam -r /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/sacCer3.bed`

Output is an xls file and a summary txt file (mean & median values across all genes in sample).

Visualise TIN in boxplots in [Rstudio](https://github.com/friedue/course_RNA-seq2015/blob/master/03_mRIN.R) using ggplot

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

__Summarising results of different QC tools with MultiQC__

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
eyJoaXN0b3J5IjpbMTA0MDYxNzU2NywtMTY4ODMwNjA2NF19
-->