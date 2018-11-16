> # Sequence Alignment

The purpose of **ASSEMBLY** is to merge reads into larger contiguous sequences (contigs) based on their sequence similarity to each other. 
The purpose of **ALIGNMENT** is to compare the contig to a previously assembled reference genome (or transcriptome) sequence.

> Mapping
> -   A mapping is a region where a read sequence is placed.
> -   A mapping is regarded to be correct if it overlaps the true region.
> 
> Alignment
> -   An alignment is the detailed placement of each base in a read.
> -   An alignment is regarded to be correct if each base is placed correctly.

For each base alignment the options are:
	- `-` a space
	- `|` match
	- `.` a mismatch

Tools used to be separated into aligners vs mappers. However these have become combined over time. However a distinction still exists. For example, studies examining SNPs and variations in a genome would be primarily alignment-oriented. However, studies focusing on RNA-Seq would be essentially mapping-oriented. The reason for this difference is that in RNA-Seq you would need to know where the measurements come from, but you would be less concerned about their proper alignment.

# Alignment Strategies

1. align reads to **genome** index to identify novel splice events. Most common approach.
2. align reads to **transcriptome** index (required transcripts to be known and annotated in the reference) - i.e. pseudoalignment. Use if you have short reads < 50bp. Using Kallisto or Salmon - much quicker but wont identify new transcripts.
4. **de novo** assembly: assemble transcripts from overlapping tags. Useful if no reference genome exists for species studied. Also if looking for complex polymorphisms that would be missed by aligning to reference. 

![enter image description here](https://lh3.googleusercontent.com/K400ZHmBCmhNY475bKN4PGdSpxK0lbqTNGBWHkWzh5DmcCuUKoDGbnuZDh6S_C_UEjPkcvTkjXIY0w "3 RNA seq mapping strategies")

Aligning to Genome vs de-novo assembly:

![enter image description here](https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/resize/user/18/Figure19-700x527.png)

# Which Aligner Should I Use?
~ 2005 high throughput short-read sequencers changed the face of sequencing which became much cheaper and accessible. Previously sequencing was laborious and expensive and the focus was on producing accurate long reads (1000bp). Sequencing reads longer improves alignment.  Sequencers now produce millions of short reads (50-300bp). Aligners have thus changed to adapt to short reads - rapidly select the best mapping at the expense of not investigating all potential alternative alignments.
- Short read aligners are incredible! They match > 10,000 sequences per second against 3 billion bases of the human genome.
- There is large variation in results between different aligners. A tool that prioritises finding exact matches will do so at the expense of missing locations and overlaps and vice versa.
- Limitations:
	- finds alignments that are reasonably similar but not exact (algorithm will not search beyond a defined matching threshold)
	- cannot handle very long reads or very short reads (<30bp) (become inefficient)

Aligners aim to:
1. speed up process
2. reduce amount of memory used
3. maximise quality of the alignment

Note the RNA aligners in red below:
![enter image description here](https://www.ebi.ac.uk/~nf/hts_mappers/mappers_timeline.jpeg)

The main RNA-seq aligner milestones are:
1. TopHat > TopHat2. Slowish 1000mins. Low memory 4GB. Everyone used it from 2009
2. STAR. Faster 24mins. Used more memory 28GB (need a server)
3. HISAT. Faster 20mins & reduce memory 4GB

|Aligner |Speed  |Memory |Year|
|--|--|--|--|
| TopHat > TopHat2 |1000mins|4GB |2009|
|STAR|24mins| 28GB |2012|
|HISAT > HISAT2 |20mins|4GB |2016|

## Splice-aware Aligners

If you are aligning to a genome then you need a splice aware aligner:
- reads span large introns
- reads represent mature mRNA have introns removed
- but aligning back to the reference genome needs to recognise intron sites
- unless reads a short <50bp then you need a splice aware aligner

![enter image description here](https://media.nature.com/lw926/nature-assets/nmeth/journal/v12/n4/images/nmeth.3317-F1.jpg)

Options:
-   [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf): Really fast, produces counts for you too. Straight forward RNA seq for differential gene expression analysis. Efficient. Sensitivie. Identified large number of novel splice sites.
-   [HISAT2](https://www.nature.com/articles/nmeth.3317): similar algorithms as Bowtie2/Tophat2 - but performs at much faster speeds
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



### **CIGAR** (Concise idiosyncratic gapped alignment report string)
CIGAR string is an alignment format used in SAM (sequence alignment map) files. CIGAR uses letters M, I, D etc to indicate how the read aligned to the reference sequence at that specific locus. In  extended CIGAR the symbol `X` is used for mismatches.
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







# Reference Genomes
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

## Reference Genome File Formats

### Reference Sequence
Reference sequences are **FASTA files.** Reference sequences are long strings of ATCGN letters.  

### Reference Annotation
Annotation files store start sites, exon, introns. They contain spacial coordinate information describing genomic loci (2D position) according to their pattern (chromosome, start, end) as intergers. These integers represent the left & right positions of an interval on the genome.

Both formats contain information:
- Span: start & end coordinates
- Hierarchy: feature belongs to another feature e.g. exon to a transcript
- Value: at a apsecific coordinate
- Annotation: feature X is also characterised by Y & Z

2 types of genome formats exist:
- [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) originally designed for UCSC genome browser
- [GFF](http://mblab.wustl.edu/GTF22.html) (also its variants GFF2, GTF & GFF3)

**BED Format**
BED is 0 based, non-inclusive on the right. This means that the interval [1, 5) contains 1,2,3,4 and coordinate 1 is the second coordinate of the genome (0 is the first).
A single line record stores all the information on the block structure of a record.

3 compulsory fields: chromosome & start & end.
9 optional fields: name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        	field number can vary from 3 - 12. Must be consistent within a file.

**GFF = General Feature Format**
GFF is 1 based, inclusive on both ends. This means that [1,5] contains 1,2,3,4,5 and that coordinate 1 is the first coordinate of the genome.
The relationship is built from data distributed over multiple lines. 

9 fields separated by TAB. 
* Reference sequence: coordinate system of annotation eg Chr1
* Source: annotation software
* Method: annotation type (eg gene)   	[GFF3 = Type: term from lite Sequence Ontology SOFA or accession number]
* Start Position: 1-based integer, <= stop position
* Stop Position: 0-length features (insertion sites)
* Score: sequence identity
* Strand: + = forward. - = reverse. “.” = no stranded
* Phase: codon phase for annotations linked to proteins. 0; 1; or 2 = indicates frame (bases to be removed from beginning of feature to reach first base of next codon)
* Group: class and ID of an annotation which is the parent of the current one.  	[GFF3 = Attributes: TAG=VALUE pairs - ID, name, alias, parent, target, gap]
* Source
* Type

Always use GFF based formats (BED is being gradually phased out).

# Download Reference Sequence & Annotation files
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


# Alignment Workflow:

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

        echo "Running STAR on $SAMPLE"
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

# Sequence Alignment Map (SAM/BAM) file

SAM/BAM files store read alignments.

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
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE4Mzc0NzIzMTksLTM3NzM0MzYxOCw5OT
g5ODg2NTYsLTE0NzA5Mjg4OTYsLTQ4Njg4NDg0NCwtMTQ3ODU2
MDQ5NiwtMTU4NjQxMzgyNiw2MzAyNDc5MDUsNjU3NTQyMjE4XX
0=
-->