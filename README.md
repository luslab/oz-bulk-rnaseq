> # RNA sequence protocol

This repository contains a protocol for the analysis of RNA-seq data. Based around the [RNA seq worksheet](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf).

__Requirements__
TBC

## Installing Bioconductor packages in R

[Install Bioconductor](https://www.bioconductor.org/install/)
[Source]("https://bioconductor.org/biocLite.R")
`biocLite (“package name“)`
For example, for erccdashboard (for artificial spike in quantification) 
`source ("https://bioconductor.org/biocLite.R")` `biocLite("erccdashboard")`

## RNA seq Workflow
1. Images: `base-calling` `demultiplexing`
2. Raw Reads `fastq` `mapping` `aligning` `STAR`
3. Aligned Reads `counting`
4. Read count table `normalise DESeq2` `edgeR`
5. List of fold changes & stats `filter`
6. Downstream analyses on differentially expressed genes

__Sequencing__
RNA extraction --> Library preparation (mRNA enrichment) --> Sequence

__Bioinformatics__ 
Process sequence reads (align) —> estimate individual gene expression levels —> normalise —> identify differentially expressed (DE) genes
 
 ## RNA Sequencing
 __RNA extraction__
* silica gel based membranes or liquid-liquid extractions with acidic phenol chloroform
* remove DNA and proteins. Improve with DNase
* quality control: Aligent bioanalyser creates an RNA integrity number (RIN) is objective way of assessing RNA quality. 10 = intact; 1 = degraded. Looks for densitometry spike at 28S and 18S rRNA bands - ratio of 28S/18S = RIN.
 
__Library Preparation__
* cDNA fragments 150-300bp —> hybridisation to flowcell (50-150 bp)
* small transcripts <150bp is lost in standard RNA-seq preparation
* mRNA enrichment: remove rRNA and tRNA by selecting polyA tails using oligodT beads OR removing rRNA with complementary sequences (ribo-minus approach —> retains unspliced RNAs).
 
__Stand sequencing__
* distinguish overlapping transcripts —> identify anti-sense transcripts by preserving which strand a fragment came from
* usually use deoxy-UTP in synthesising the 2nd cDNA strand
* hybridise DNA fragments to flowcell via adapters —> clonal amplify fragments forming clusters of dsDNA = improve signal of each fragment
* Illumina seuqencing protocols: covers 50 - 100bp of each fragment
* fragment ends is based on labelled dNTPs with reversible terminator elements —> incorporated + excited by a laser —> enables optical identification of bases
* coverage = number of reads sequences in relation to genome size (how many times each base of the genome is referenced) - for RNA-seq the size of the transcriptome is not accurately known
* Lander-waterman equation: coverage = (read length + read number)/haploid genome length
* every base should be covered more than once to identify sequence errors
* coverage is not uniform: euchromatin is overrepresented, GC rich regions are favoured by PCR
* for RNA-seq use the least abundant RNA species of interest to determine the number of required reads (= sequence depth)

Estimating the sequence depth (capture enough fragments of the least expressed genes) = ENCODE guidelines
* experiment type and biological question
* transcriptome size
* error rate of the sequencing platform

Deeper sequencing required for:
* low expressed genes
* identify small changes between conditions
* quantify alternative splicing
* detect chimeric transcripts; novel transcripts; start and end sites

Differential gene expression analysis: prioritise increasing the number of biological replicates rather than sequencing depth
* Single read vs. paired end read:
* single read = determines the DNA sequence of just one end of each DNA fragment
* paired end = yields both ends of each DNA fragment. more expensive but increase mappability for repetitive regions —> easier to identify structural variations + indwells
* for detecting de novo transcriptome assembly in humans need 100-200 x10^6 paired end reads.

## Bioinformatics
__Experimental Design__
 
Variability in results
* need replicates to capture breadth of isolate noise
* Technical replicates: repeat library preparations from the same RNA sample —> avoid batch effects, lane effects. Should multiplex same sample over different lanes of same flowcell
* Biological replicates: parallel measurements on different samples = RNA from independent cells/tissues. Most RNA seq have 3 biological replicates but ideally need 6 per condition
 
Artificial RNA spike-ins
* used to accurately quantify absolute transcript concentration. RNA of known quantities used for calibration eg ERCC. R package = erccdashboard. Different spike-in controls are needed for each RNA type.
* dont use spike ins to normalise between different samples (they dont account for differences in amount of starting material).
 
Blocking and randomise:
* Randomly choose which samples to treat and sample
* Block samples into groups based on known sources of variation (sex, weight, cell cycle status) - subexperiments in each block increases sensitivity.
 
### Raw Data (Sequencing Reads)
SRA = sequencing read archive = main repository for nucleic acid sequences (includes USA NCBI + European Bioinformatics Institute + DNA Databank of Japan)

* FASTQ files are the formate in which we store sequencing reads. There are lots of other formats but FASTQ is the most common.
* FastQ files end in SRR and SRX
*FASTQ files bundle the sequence of each single read with the quality score.
 
__Download FASTQ files__
[ENA](https://www.ebi.ac.uk/ena) OR  [SRA](https://www.ncbi.nlm.nih.gov/sra)
Search accession number (indicated in published paper)
DOWNLOAD
click link in column (Fastq files ftp) & save OR copy link address of Fastq files column —> command line & move target directory: `$ wget <link copied from the ENA website >`
Alternatively if there are many samples then download summary (right click on TEXT) & copy link location: `$ wget -O samples_at_ENA .txt "<LINK copied>"` NB the quotation marks are crucial 
Change directory `cd` to where you will store data & use 11th column of TEXT file (Fastq file top) to feed the URLs of different samples `$ cut -f11 samples_at_ENA . txt | xargs wget`
 
__Fastq-dump NCBI tool to convert fastq.sra files__
fastq.gz  = compressed version of fast file (needs unzipping)

[NCBI Format Guide](https://www.ncbi.nlm.nih.gov/books/NBK242622/)

FASTQ files are uncompressed & large. They contain:
1. `@ then read ID +- informs on sequencing run`
2. `Sequenced bases`
3. `+ (read ID again or description)`
4. `Quality score of each base (ASCII-encoded)`
 
Read ID format:  
`@<machine_id>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos><read>:<is filtered>:<control number>:<index sequence>`
 
__Paired End Sequencing__
* 2 FASTQ files: 1 forward read; 1 backward read
* know the origin of each read (forward vs reverse) - encoded in read name - some analysis tools require combining the 2 files into 1.
* Need to process the read as Split Reads/files
	* Look at to explain process `fastq help`
* The forward read will usually be filename_1 and backward read is filename_2
 
### Quality Scores
* The first bioinformatic step is quality control. Use `fastqc` - see help page by typing `fastqc -h`
* Base calling = deduce the nucleotide letter code sequence from the fluorescence signal edited when incorporated into the sequence read. Imperfect.
* Phred score, Q = proportional to probability that a base call is incorrect. 10 = 1 in 10 bases are wrong (90% accuracy); 20 = 1 in 100 bases are wrong (99% accuracy). Higher Phred = higher quality
* Sanger also have a quality score
* ASCII character = represents that Phred Score. Depends on:
	* sequence technology used
	* base caller assignment used (eg Bustard, RTA0 HiSeq X). 
* Maximum score is 45.
* converting Illumina FASTQ file 1.3 (Phred + 64) to version 1.8 (Phred +33) use: 

`sed -e '4~4y/ @ABCDEFGHIJKLMNOPQRSTUVWXYZ [\\]^_abcdefghi/!"#$%& '\ ' '()*+ , -.\/0123456789:; <= >? @ABCDEFGHIJ /' originalFile.fastq`
[Base calling](https://academic.oup.com/bib/article/12/5/489/268399)
 
If the quality scores contain character 0 it is either Sanger phred+33 or Illumina 1.8+ phred+33. When they also contain the character J, it is Illumina 1.8+ phred 33, otherwise it is Sanger phred + 33. 
When the quality scores do not contain 0, it is either Solexa +64, Illumina 1.3+ Phred+64, Illumina 1.5+ Phred+64.
It is Illumina 1.3 phred + 64 when it contains A It is Illumina 1.5 phred +64 
 
### Quality Control
Main points for QC:
* FASTQC on raw sequenced reads
* RSeQC on aligned reads
* Descriptive plots - assess read counts
 
__1st quality control point__ = RAW READS = FastQ file —> FastQC program  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
`wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip`
`unzip fastqc_v0.11.7.zip`
`cd FastQC`
`chmod 755 fastqc`

Check for:
* PCR duplicates
* adapter contamination
* rRNA + tRNA reads
* unmappable reads (contaminating nucleic acids) - FASTQC doesn’t check this
* each test will either = pass; warn; or fail. Fail is expected in some cases and does not mean the sequencing needs repeated.
 
__run FastQC__
`mkdir fastqc_results` # make a folder to store the results
run FastQC (for the course it 's available in the software folder )
`/home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/FastQC/fastqc ERR458493.fastq .gz -- extract -o fastqc_results`
Look at the results
`ls fastqc_results/ERR458493_fastqc/`
`cat fastqc_results/ERR458493_fastqc/summary .txt`
 
To summarise multiple QC outputs use the [MultiQC tool](http://multiqc.info/)

run FastQC on all fastq.gz files per sample
$ for SAMPLE in WT_1 WT_2 WT_3 WT_25 # random selection of samples
`mkdir fastqc_results/${SAMPLE}`
run `multiqc` within the `fastqc_results` folder and use the folder names (WT_1 etc.) as prefixes to the sample names in the final output
 
## Read Alignment / Mapping
Identify cDNA transcripts in a sample by mapping them to the genomic origin using a reference genome
Aim is to map millions of reads accurately and quickly
Limitations are sequencing errors; genomic variation; repetitive elements
Main Challenge of RNA seq = spliced alignment of exon-exon spanning reads (eg Read 2 in image); multiple different transcripts (isoforms) from same gene

__RNA seq Programmes (STAR, TopHat, GSNAP)__
        	1. align reads to transcriptome (required transcripts to be known and annotated)
        	2. identify novel splice events (using reads that cant be aligned to reference transcriptome)

* False positives: lowly expressed isoforms are excluded by algorithms —> bias towards identifying strongly expressed genes
* Mapping ambiguity = many reads overlap with more than one isoform
* Sequencing reads longer improves alignment
* Alignment-free transcript quanitification = ignore location of reads in a transcript; compare k-mers of the reads in hash tables of transcriptome and note matches.
 
### Reference Genomes
As a rule for human data use GenCODE. For other species use Ensemble.
ENCODE, iGenomes, NCBI, UCSC, Mouse Genome Project, Berkeley Drosphilia Project

Reference sequences = FASTA files. 
Compress with `gzip` command or `faToTwoBit` (into .2bit files)\

Reference sequences are long strings of ATCGN letters. 
File formats store start sites, exon, introns. One line per genomic feature.
 
__File Formats__
GFF = General Feature Format
9 fields separated by TAB. No white space - all fields except last in each feature must contain value (missing values = ‘.’)
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
 
GTF = Gene Transfer Format (aka GFF 2.5). More strict than GFF. Same as GFF2 but 9th field expanded into attributes (like GFF3). http://mblab.wustl.edu/GTF2.html

__Download GTF & FASTA files__
UCSC https://genome.ucsc.edu/          	https://genome.ucsc.edu/cgi-bin/hgTables
ENSEMBL http://www.ensembl.org/index.html
RefSeq
GenCODE

* UCSC and Ensembl use different naming conventions (which impacts on analyses) - try to stick to one.
 
Convert 2bit format —> FASTA format: 
`twobittofa file_name.2bit file_name.fa`

_from ENSEMBL:_
http://www.ensembl.org/info/data/ftp/index.html
Search for species Saccharomyces cerevesiae
Click on Gene sets GTF link & DNA FASTA link
GTF: Right click on Saccharomyces)cerevisiae.R64-1-1.92.gtf.gz → copy link address
FASTA: Right click on DNA top level file.
In command line (in appropriate Folder) `wget [paste link address]`
Unzip file `gunzip file_name`

_From UCSC:_
Download a GTF file of yeast transcripts from the UCSC Genome Table Browser https://genome.ucsc.edu/cgi-bin/hgTables
Move the downloaded gtf file to the appropriate folder
**![](https://lh4.googleusercontent.com/piQkvTkiSIYCY9m-gATKN8CTmWGFPVZaP7KItC44zJP_oztaNMxjf9O33hljoHvARnSAqaXP1lz5pUo8_7X49xlHKXtX5hUyU-vAfehxNnXAVQ3mh152qUNwlywheUpx5P2GUa4Y)**
N.B. GTF files downloaded from UCSC table have same entries for gene ID and transcript ID → creates problem with analysing different exon isoforms (same gene ID but different transcriptt ID)

__Install genePredToGft__
`wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
Linux version - download to bin folder
`bash Miniconda3-latest-Linux-x86_64.sh  https://conda.io/docs/user-guide/install/linux.html`
`wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/`
`source .bash_profile`
`conda install -c bioconda ucsc-genepredtogtf `
remove first column & first line: `cut -f 2- file_name.txt | sed ‘1d”`
`genePredToGtf file file_name file_name.gtf`

__BED Format simplest annotation store__
3 compulsory fields: chromosome & start & end.
9 optional fields: name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        	field number can vary from 3 - 12. Must be consistent within a file.
indicates region with 0-based start and 1-based end position (GFF & GTF are 1-based in both start and end)
Aligning Reads
 
### 1. Choose alignment tool
 
* Multiple alignment programmes available, each specialising in detecting different factors eg structural variants; fusion transcripts
* Straight forward RNA seq for differential gene expression analysis = use STAR [STAR manual PDF](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
	* Efficient
	* Sensitive
	* But large number of novel splice sites (caution)
* TopHat = popular aligner – wrapper around the genomic aligner Bowtie
* The alignment tool has relatively little impact on the downstream analyses (vs. annotation, quantification, differential expression tools)
 
### 2. Generate  Genome Index
 
Input Files = Reference Genome & Annotation File

![RNA-seq Flowchart - Module 1](https://github.com/griffithlab/rnaseq_tutorial/wiki/Images/RNA-seq_Flowchart2.png)

Genome sequence
Suffix Arrays
Chromosome names & lengths
Splice junction coordinates
Gene information
Create directory to store index in: `mkdir STARindex`
Module load STAR `ml STAR`  

STAR command line has the following format:
`STAR --option1-name option1-value(s)--option2-name option2-value(s) ...`
If an option can accept multiple values, they are separated by spaces, and in a few cases - by commas.

Send cmd to generate genome as `batch job` to cluster:
`sbatch -N 1 -c 8 --mem 32G --wrap="STAR --runMode genomeGenerate --genomeDir /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/STARindex --genomeFastaFiles /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --sjdbGTFfile /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/Saccharomyces_cerevisiae.R64-1-1.92.gtf --sjdbOverhang 49 --runThreadN 8"`

`N` is the number of nodes (usually leave at 1)
`c` is the number of cores (i.e. threads - so you could change the STAR command to `runThreadN 8` with this example)
`--mem` is the amount of memory you want

no space between -- and the parameter
quotation marks: use ` and not ' or "

no `<>` needed
must use Ensembl FASTA file with Ensemble GTF file. You cannot mix Ensembl with UCSC without modifying first and this isn't advised as they have different scaffolds.

Check it is running using `myq`
If it isn’t seen there then `ls` the folder → look for output file named “slurm…”
Open slurm file: `more slurm…` which will explain outcome of file eg FATAL INPUT PARAMETER ERROR
 
Alternative approach (assign runSTAR & REF_DIR variables):
`${runSTAR } --runMode genomeGenerate \ --genomeDir ~/STARindex\` # index will be stored there
 `--genomeFastaFiles ${REF_DIR}/sacCer3 .fa\` # reference genome sequence
`--sjdbGTFfile ${REF_DIR}/sacCer3.gtf\` # annotation file
`--sjdbOverhang 49` # should be read length minus 1 ; length of the genomic sequence around the annotated junction to be used for the splice junctions database
`--runThreadN 1\` # can be used to define more processors
 
### 3. Align each FASTQ file

* Need to align each FASTQ file
* Sample distributed over n = X flow cell lanes → X fastq files per sample
* STAR merges the X files if multiple file names are indicated (other align tools dont)
* Separate the file names with a comma (no spaces)
* Create directory to store STAR output `mkdir alignment_STAR`
* List fast.qz files separated by commas and remove white spaces:

INPUT=`ls -m /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/SNF2_rep1/*.fastq.gz| sed 's/ //g' | echo $INPUT | sed 's/ //g'`

`ls -m` list, fill width with a comma separated list of entries all the `fast.gz` files = compressed filed.
`sed` remove a space from each file name
`echo` displays a line of text - in this case it lists all the file names  - in this case the pipe to echo is to remove new spaces that are created between multiple lines.
no space after FILES - with space after it thinks FILES is command.
`sed` = stream editor - modify each line of a file by replacing specified parts of the line. Makes basic text changes to a file - `s/input/output/g`

Execute STAR in runMode default `alignReads`
For WT_rep1 folder (already uncompressed):
`sbatch -N 1 -c 8 --mem 32G --wrap="STAR --genomeDir /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/STARindex --readFilesIn /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/WT_rep1/$INPUT --outFileNamePrefix /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/WT_1_ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 1"`

For SNF2_rep1 folder (compressed):
`sbatch -N 1 -c 8 --mem32G --wrap="STAR --genomeDir /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/STARindex --readFilesIn /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/SNF2_rep1 --readFilesCommand gunzip -c  --outFileNamePrefix /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/SNF2_1_ --outFilterMutlimapNmax 1 \ --outReadsUnmapped Fastx \ --outSAMtype BAM SortedByCoordinate \ --twopassMode Basic \--runThreadN 1`

Check the summary of the output:
`cat file_name_Log.final.out`
![enter image description here](https://user-images.githubusercontent.com/33317454/37438378-095cb442-282d-11e8-95fd-bfecefae5b75.png)

By allocating all file names to the `$INPUT` term it means all the FASTQ files are read together but it is best to do each separately and use the `snakemate` command - this makes it easier to see if there is an error in an individual sequencing file.

The STAR PDF manual has all the explanations on how to write the STAR command and fine tune parameters e.g.:
* multi-mapped reads
* optimise for small genomes
* define min & max intron sizes 
* handle genomes with >5000 scaffolds
* detect chimeric & circular transcripts

ENCODE options:
'outFilterMultimapNmax 1' max number of multiple alignments allowed for a read - if exceeded then read is considered unmapped i.e. when 1 is used this only identifies unique reads and removes multimapped reads. This is generally accepted.

#`STAR` will perform mapping , then extract novel junctions which will be inserted into the genome index which will then be used to re -map all reads
#`runThreadN` can be increased if sufficient computational power is available

### 4. BAM file indexing

Create BAM.BAI file with every BAM file to quickly access the BAM files without having to load them to memory
Install samtools

Run samtools index cmd for each BAM file once mapping is complete:
`ml SAMtools`
`samtools index /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam`

### 5. Store Aligned Reads: SAM/BAM files

__SAM Files__

SAM file = Sequence Alignment Map - generic nucleotide alignement format describing alignment of sequenced reads to a reference. [More Details](https://github.com/samtools/hts-specs) here.
* Contain short header & long alignment sections
* Each row represents a single read alignment
	* starts with @ then abbreviation: SQ = sequence directory listing chromosomes names (SN) and lengths (LN)
* Each read has 11 mandatory entries (black font) & optional fields (grey font)

![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/bam_structure.png)

_SAM header section_

* Begin with @, followed by `tag:value pairs`
	* `tag` is two-letter string defining the value e.g. `@SQ` = names & lengths of reference sequences
* 1 line per chromosome
* To retrieve only the SAM header `samtools view -H`
* To retrieve both the header & alignment sections `samtools view -h`
* The default `samtools view` will not show the header section

`samtools view -H WT_1_Aligned.sortedByCoord.out.bam`

_SAM alignment section_

* Each line coresponds to one sequenced read
* 11 mandatory fields in specified order:

`<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>`

* If info is unavailable then value can be 0 or * 
* Optional fields follow the mandatory fields

![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/sam_fields.png)

FLAG field:
* stores info on the respective read alignment in one single decimal number
* decimal is the sum of all the answers to Yes/No questions:
![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/sam_flag.png)
To convert the FLAG integer into plain english [click here](https://broadinstitute.github.io/picard/explain-flags.html).

`CIGAR` = Concise idiosyncratic gapped alignment report string
* Indicates which operations were necessary to map the read to the reference sequence at that specific locus
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

_Optional fields_
Following the 11 mandatory fields, the optional fields are presented as key-value pairs in the format of  `<TAG>:<TYPE>:<VALUE>`, where  `TYPE`  is one of:
-   `A`  - Character
-   `i`  - Integer
-   `f`  - Float number
-   `Z`  - String
-   `H`  - Hex string

Reads within the same SAM file may have different numbers of optional fields, depending on the program that generated the SAM file. 

Commonly used optional tags include:
-   `AS:i`  - Alignment score
-   `BC:Z`  - Barcode sequence
-   `HI:i`  - Match is i-th hit to the read
-   `NH:i`  - Number of reported alignments for the query sequence
-   `NM:i`  - Edit distance of the query to the reference
-   `MD:Z`  - String that contains the exact positions of mismatches (should complement the CIGAR string)
-   `RG:Z`  - Read group (should match the entry after ID if @RG is present in the header.

E.g. we can use the NM:i:0 tag to select only those reads which map perfectly to the reference (i.e., have no mismatches). Tags that begin with  `X`,  `Y`, and  `Z`  are reserved for particularly free usage and will never be part of the official SAM file format specifications.  `XS`, for example, is used by TopHat to encode the strand information (e.g.,  `XS:A:+`) while Bowtie2 and BWA use  `XS:i:`  for reads with multiple alignments to store the alignment score for the next-best-scoring alignment (e.g.,  `XS:i:30`).

BAM file = Binary Alignment Map - human readable TAB-delimited compressed. Bigger than gzipped SAM files as they are optimised for rapid access (not just size reduction).
Position sorted BAM files = indexed so all reads aligning to a locus can be retreived without loading the entire file to memory.

__Convert a BAM file to SAM file__
`samtools view -h WT_1_Aligned.sortedByCoord.out.bam > WT_1_Aligned.sortedByCoord.out.sam`

__Read Group__
Key feature of SAM/BAM format is ability to label individual reads with readgroup tags. This allows pooling results of multiple experiments into a single BAM dataset to simplify downstream logistics into 1 dataset.
Downstream analysis tools eg `variant callers` recognise readgroup data and output results.

[BAM readgroups GATK website](https://gatkforums.broadinstitute.org/gatk/discussion/1317/collected-faqs-about-bam-files)

Most important readgroup tags: `ID`, `SM`, `LB`, `PL`
![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/rg.png)

## Manipulating SAM/BAM files
There are 4 major toolsets for processing SAM/BAM files:

- [SAMTools](http://www.htslib.org/) - interact with high throughput sequencing data, manipulate alignments in SAM/BAM format, sort, merge, index, align in per-position format
- [Picard](https://broadinstitute.github.io/picard/) - Java tools for manipulating high throughput sequencing data
- [DeepTools](https://deeptools.readthedocs.io/en/develop/) - visualise, quality control, normalise data from deep-sequencing DNA seq
- [BAMtools](https://github.com/pezmaster31/bamtools/wiki/Tutorial_Toolkit_BamTools-1.0.pdf) = read, write & manipulate BAM genome alingment files

__Processing:__
1. Filter using BAM Tools
	-  mapping quality: remove poor  alignments - eg remove all alignments below Phred scale quality of 20
	- keep only those which are "properly paired" ie forward is looking at reverse (for paired reads)
	- reference chromosome: e.g. only keep mitochondrial genome alingments
2. Remove duplicates with Picard
3. Clean up with CleanSam Picard tool
	- fixes alignments that hang off ends of ref sequence
	- sets MAPQ to 0 if read is unmapped

To peak into a SAM or BAM file: `samtools view FILENAME.bam | head`
Convert a BAM file into a human readable SAM file (including the header): `samtools view -h FILENAME.bam > FILENAME.sam`
Compress a SAM file into BAM format (-Sb = -S -b)" `samtools view -Sb FILENAME.sam > FILENAME.bam`
Generate an index for a BAM file (needed for downstream tools): `samtools index FILENAME.bam`

SAMTools help page = `samtools --help`
Usage:   `samtools <command> [options]`

5 key commands:
1. Indexing
2. Editing
3. File operations (aligning, converting, merging)
4. Statistics
5. Viewing

For each command there are multiple options ` samtools COMMAND -X` 
Create BAM with only **unmapped reads**: `samtools view -h -b -f4 FILENAME.bam > unmapped_reads.bam` 
Create BAM with only **mapped reads**`samtools view -hb -F 4 FILENAME.bam > mapped_reads.bam` 
Create BAM with **mapping quality >= 20**`samtools view -h -b -q 20 FILENAME.bam > high_mapq_reads.bam` 
Create BAM with **uniquely aligned reads** (STAR gives uniquely aligned reads a mapping quality of 255 so you can use samtools to pull all reads with mapping quality = 255 only (using samtools command, option -q, = 255) `samtools view -h -q 255 FILENAME.bam > uniquely_aligned_reads.bam`

- `-h` is used to print the header (always needed).
- sam files are human readable, bam are compressed. BAM are much smaller
	- `-b`will produce a BAM file. `-s` will produce a SAM file.
	- for SAM files you can run other commands on them eg head FILENAME.sam whereas BAM files need to be run through samtools i.e. `samtools view FILENAME.bam | cut -f 2 | head`

Create BAM with only **reads aligned to reverse strand**:
- First, sort the BAM file (A-Z; 0-9) using `sort`
- Second, count the adjacent lines which are identical using `samtools view FILENAME.bam | cut -f 2 | sort | uniq -c`
		- Third, create file with the specific feature e.g. reverse reads (FLAG = 16 in column 2): `samtools view -h -f 16 FILENAME.bam > reverse_reads.bam`
		
Alternatively used awk command using SAM file:
- `cut -f 2 FILENAME.bam | sort | uniq -c`
- `awk '{OFS="\t"} $1 ~ /^@/ || $2==0 {print $0}' FILENAME.sam | cut -f 2| sort | uniq -c`

where:
`cut -f 2` selects out only the field list 2 (2nd column)
`{OFS="\t"}` converts spaces to tabs as awk would ordinarily use spaces but the file is separated by tabs
`$1 ~ /^@/` represents the header line where $1 = column 1, ~ means looks like, ^@ means start with @
`||` means OR
`$2==0` means column 2 equals 0 i.e. FLAG (column 2) is forward strand (forward strand = 0, reverse strand = 16) - decode all SAM file numbers with https://broadinstitute.github.io/picard/explain-flags.html
`{print $0}` means print everything ($0 means everything)

- To **filter** an alignment file using the **optional tags** you have to use other tools eg `grep` to look for exact matches.
- Optional SAM/BAM fields depend on the alignment program used - before filtering make sure you know how the aligned generated the value.

Create SAM file with **reads of insert sizes > 1000bp** use the CIGAR string (column 6 in SAM file):
- first convert BAM to SAM file ` samtools view -h FILENAME.bam > FILENAME.sam`
- use `grep` to exclude (using `-v`) lines with >3 digits (using `[0-9][0-9][0-9][0-9]`) followed by `N` (N means mismatch i.e. skipped bases, whereas M = match) `egrep -v "[0-9][0-9][0-9][0-9]N" FILENAME.sam > smallinsert_reads.sam`
- Alternatively use `awk` to focus on column 6 (`$6` in CIGAR string) and exclude lines with 3 digits (using `![0-9][0-9][0-9][0-9]`) then printing everything `{print $0}` and creating new file:  `awk '!($6 ~ /[0-9][0-9][0-9][0-9]N/) {print $0}' FILENAME.sam > smallinsert_reads.sam` 

Create SAM file with **intron spanning reads**:  
- use `grep` to select lines with a number of digits (using `[0-9]+`) then `M` (i.e. matches) then any number of digits again, then `N` (i.e. mismatches) then any number of digits and then M again at the end: `egrep "(^@|[0-9]+M[0-9]+N[0-9]+M)" FILENAME.sam > intron-spanning_reads.sam`
- Alternatively use awk to focus on column 6 (CIGAR string) and select the header `$1 ~ /^@/` and the 6th column with any number of digits followed by M followed by digits then N then digits the M: `awk '$1 ~ /^@/ || $6 ~ /[0-9]+M[0-9]+N[0-9]+M/ {print $0}' FILENAME.sam > intron-spanning_reads.sam`

## Quality control of aligned reads
**STAR output files:**
 - Aligned.sortedByCoord.out.bam - the loci of each read & sequence
 - Log.final.out - summary of alignment statistics
 - Log.out - commands, parameters, files used
 - Log.progress.out - elapsed time
 - SJ.out.tab - loci where splice junctions were detected & read number overlapping them
 - Unmapped.out.mate1 - fastq file with unmapped reads

More information on these is in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) chapter 4 page 10.
 
 The main STAR output file for downstream analyses are
 **Aligned.sortedByCoord.out.bam** & **Log.final.out**.

**out.mate1 files** are fastx files that are uncompressed and can be made smaller using gzip

## Quality control of aligned reads

After aligning and before performing downstream analyses check for:
1. Excessive amounts of reads not aligned
2. Obvious biases in the read distributions
3. Similarity between replicate samples

__Alignment Assessments__
 
 **1. Check that mapping rate of RNA-seq reads is > 70%**
 
 Check the aligner's output: `cat FILENAME_Log.final.out`
 - most important number = **uniquely mapped reads**
 - if using >2 BAM files then visualise alignment rate for each e.g. using MultiQC in R studio: see https://github.com/friedue/course_RNA-seq2015/blob/master/01_Alignment_visualizeSTARresults.pdf
Mount files onto laptop: Right click on Finder --> Connect to server --> Connect to Luscombe Lab
- `infiles <- list.files(path="/Volumes/lab-luscomben/working/oliver/projects/rna_seq_worksheet/alignment_STAR", pattern = "WT_1_Log.final.out", full.names = TRUE)`
- `align.results <- lapply(infiles, function(x) read.table(x, sep="|", strip.white = TRUE, stringsAsFactors = FALSE, skip = 3, fill = TRUE, header = FALSE))`
- `typeof(align.results)`
- `head(align.results[[1]])`
- > `align.results <- lapply(align.results, function(x) transform(x,V2 = as.numeric(gsub("%", "", x$V2) )))`


 **2. Calculate number of alignments in each BAM file**
- easiest to do using a line count ` samtools view Aligned.sortedByCoord.out.bam | wc -l`
- unmapped reads in the BAM file & also multiple instances of the same read mapped to different locations will also be counted (latter only if multi-mapped reads were kept) so run specific tools to indicate FLAG values too.
- `samtools flagstat` assesses the FLAG field and prints a summary report: `samtools flagstat Aligned.sortedByCoord.bam`


## Visualising Data


__Visualise STAR alignment information in R studio using ggplot2__
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

__Visualise the output of  RSeQC quality controls__
Basic alignment stats: `bam_stat.py -i WT_1_Aligned.sortedByCoord.out.bam`

To add results of samtools flagstat & RSeQC to a MultiQC report capture the output as a txt file.
`bam_stat.py -i WT_1_Aligned.sortedByCoord.out.bam > bam_stat_WT_1.txt`
`samtools flagstat WT_1_Aligned.sortedByCoord.out.bam > flagstat_WT_1.txt`

To visualise the output of mulple RSeQC reads download the relevant txt files and follow this [R script](https://github.com/friedue/course_RNA-seq2015/blob/master/02_Alignment_QC_visualizeReadDistributionsAsBarChart.R).

__Visualise Aligned Reads__
Check results visually to ensure reads align to expected regions without excess mismatches.

Genome Browsers:
- Broad institute [Integrative Genomics Viewer IGV](http://software.broadinstitute.org/software/igv/book/export/html/6)  
- Ensembl
- UCSC

![IGV image](http://software.broadinstitute.org/software/igv/sites/cancerinformatics.org.igv/files/images/igv_desktop_callouts.jpg)

IGV = view reads in visual format
http://software.broadinstitute.org/software/igv/
shows the transcripts amount related to an annotated genome
there is often mismatch between different annotations eg ref-seq and gencode: choosing between the two is controversial
top row = chromosome. red bar is location. blue lines mid-section refer to transcripts binding with more = higher peak. bottom section = reference genomes.

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
