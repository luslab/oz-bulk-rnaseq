<<<<<<< HEAD
=======
# RNA_seq_protocol_OZ
# RNA sequence protocol

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
ENCODE, Mouse Genome Project, Berkeley Drosphilia Project

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
 
1. Choose alignment tool
 
* Multiple alignment programmes available, each specialising in detecting different factors eg structural variants; fusion transcripts
* Straight forward RNA seq for differential gene expression analysis = use STAR [STAR manual PDF](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
	* Efficient
	* Sensitive
	* But large number of novel splice sites (caution)
* TopHat = popular aligner – wrapper around the genomic aligner Bowtie
* The alignment tool has relatively little impact on the downstream analyses (vs. annotation, quantification, differential expression tools)
 
2. Generate Genome Index
 
Input Files = Reference Genome & Annotation File
Genome sequence
Suffix Arrays
Chromosome names & lengths
Splice junction coordinates
Gene information
Create directory to store index in: mkdir STARindex
Module load STAR `ml STAR`  
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
 
3. Align each FASTQ file

Sample distributed over n = X flow cell lanes → X fastq files per sample
STAR merges the X files if multiple file names are indicated (other align tools dont)
Separate the file names with a comma (no spaces)
Create directory to store STAR output `mkdir alignment_STAR`
List fast.qz files separated by commas (no spaces):
`$FILES =ls -m /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/WT_rep1/*.fastq.gz| sed 's/ //g'`
`$FILES = echo $FILES | sed 's/ //g'`

`sed` = stream editor - modify each line of a file by replacing specified parts of the line. Makes basic text changes to a file


Execute STAR in `runMode “alignReads”`
`${ runSTAR } --genomeDir ${REF_DIR}/STARindex/ --readFilesIn $FILES --readFilesCommand zcat \ --outFileNamePrefix alignment_STAR/WT_1_ --outFilterMultimapNmax 1 \ --outReadsUnmapped Fastx \ --outSAMtype BAM SortedByCoordinate \ --twopassMode Basic \ --runThreadN 1`

#necessary because of gzipped fastq files
#only reads with 1 match in the reference will be returned as aligned
#will generate an extra output file with the unaligned reads
#`STAR` will perform mapping , then extract novel junctions which will be inserted into the genome index which will then be used to re -map all reads
#`runThreadN` can be increased if sufficient computational power is available

4. BAM file indexing

Create BAM.BAI file with every BAM file to quickly access the BAM files without having to load them to memory
Install samtools
Run samtools index cmd for each BAM file once mapping is complete
export PATH =~/working/oliver/bin/samtools -1.7: $PATH
samtools index alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam



## Visualising Transcripts
 
IGV = view reads in visual format
http://software.broadinstitute.org/software/igv/
⁃                shows the transcripts amount related to an annotated genome
⁃                there is often mismatch between different annotations eg ref-seq and gencode: choosing between the two is controversial
⁃                top row = chromosome. red bar is location. blue lines mid-section refer to transcripts binding with more = higher peak. bottom section = reference genomes.
>>>>>>> f197036ca7a92787d1aefd98f22855d6ddfbb312
