> # RNA sequence protocol assessing for Alternative Splicing & Polyadenylation

- This repository contains a protocol to analyse RNA-seq data, focusing on alternative splicing & polyadenylation, authored by Oliver Ziff. 
- The contents are based on multiple resources including the [RNAseq worksheet](http://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf); [Biostars handbook](https://www.biostarhandbook.com/); [Data Camp](https://www.datacamp.com/home); [Coursera](https://www.coursera.org/specializations/bioinformatics); and most importantly the experience of established experts in RNAseq analysis within [my host laboratory](https://www.luscombelab.org/crickmembersdetail). 
- The protocol utilises a combination of bash `unix` commmand line and `R` scripts.

# RNA-seq Workflow
This forms the chapters in this repository:

**Wet-lab sequencing phase:**
1. Extract & isolate RNA
2. Prepare library: break RNA into small fragments, convert to dsDNA, add sequencing adapters, PCR amplify
3. Strand Sequence the cDNA library: flow cell, base calling & quality score, replicates (technical = multiple lanes in flow cell; biological = multiple samples from each condition)
![Preparing RNA seq library](https://lh3.googleusercontent.com/RYpyReGfJbJOWjm20hzclqR6KUMkacZ6p_xaKvQs3piOTfxXdRiXUmiKAd45nHWj30cxJPVXmqTfnQ)
![enter image description here](https://lh3.googleusercontent.com/EBRN0O87F248JvjOzL_yHF1U328THjmXywtF4shxKxmzIwePgU-XR6ETv9Q0LCFP7bEcltsTXrN9hg)

**Bioinformatic phase:**
1. Experimental design: variability, spike-ins, blocking & randomise, filter out low quality reads & artifacts (adapter sequence reads)
2. Raw Reads: FATQ files download SRA, quality scores (Phred), paired vs single end sequence, FASTQC quality control
3. Align (map) reads to reference genome (FASTA, GFF, GTF): annotation file (BED), alignment program (STAR, HISAT), reference genomes (GenCODE, Ensemble), generate genome index, create & manipulate BAM/SAM files containing sequence alignment data
4. Visualise alingment data in R studio: ggplot2, IGV genome browser, sashimi plots, bias identification QoRTs, read quantification with gene based read counting 
5. Normalise between samples & Log Transform read counts: adjust each gene read counts for the total aligned reads in within each sample. Log2 scale, visually explore, variance shrinkage, 
6. Plot the data using global read count patterns: as there are 20,000 genes with multiple samples there are too many data points to plot everything. Summarise data with pairwise correlation, hierarchical clustering, PCA analysis - look for differences between samples & identify outliers to consider excluding
7. Identify differentially expressed genes between normal & mutants:  poisson distribution, exploratory plots (histograms, MA plot, heatmaps, read counts of single genes). Use `DESeq2` `edgeR`

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
 
 # Wet-lab RNA Sequencing
 __RNA extraction__
* silica gel based membranes or liquid-liquid extractions with acidic phenol chloroform
* remove DNA and proteins. Improve with DNase.
* quality control: Aligent bioanalyser creates an **RNA integrity number (RIN)** is objective way of assessing RNA quality & degradation. 10 = intact; 1 = degraded. RIN of 8 is generally accepted threshold before proceeding to RNA seq. Uses elecetophoresis and looks for densitometry spike at 28S and 18S rRNA bands - ratio of 28S/18S = RIN.
 ![enter image description here](http://tlcr.amegroups.com/article/viewFile/286/596/2055)
__Library Preparation__
* cDNA fragments 150-300bp —> hybridisation to flowcell (50-150 bp)
* small transcripts <150bp is lost in standard RNA-seq preparation
* mRNA enrichment: remove rRNA and tRNA by selecting polyA tails using oligodT beads OR removing rRNA with complementary sequences (ribo-minus approach —> retains unspliced RNAs).
 
__Stand sequencing__
* distinguish overlapping transcripts —> identify anti-sense transcripts by preserving which strand a fragment came from
* usually use deoxy-UTP in synthesising the 2nd cDNA strand
* hybridise DNA fragments to flowcell via adapters —> clonal amplify fragments forming clusters of dsDNA = improve signal of each fragment
* Illumina seuqencing protocols: covers 50 - 100bp of each fragment
* fragment ends is based on labelled dNTPs with reversible terminator elements —> incorporated & excited by a laser —> enables optical identification of bases
* coverage = number of reads sequences in relation to genome size (how many times each base of the genome is referenced) - for RNA-seq the size of the transcriptome is not accurately known
* Lander-waterman equation: coverage = (read length + read number)/haploid genome length
* every base should be covered more than once to identify sequence errors
* coverage is not uniform: euchromatin is overrepresented, GC rich regions are favoured by PCR
* for RNA-seq use the least abundant RNA species of interest to determine the number of required reads (= sequence depth)

Estimate the sequence depth (aim is to capture enough fragments of the least expressed genes)
-  Recommendations from [ENCODE guidelines](https://www.encodeproject.org/about/experiment-guidelines/)
- experiment type and biological question
- transcriptome size
- error rate of the sequencing platform

Deeper sequencing of RNA is required to:
* identify low expressed genes
* identify small changes between conditions
* quantify alternative splicing (intron retention & exon skipping)
* detect chimeric transcripts; novel transcripts; start and end sites

Prioritise increasing the number of biological replicates rather than the sequencing depth

Illumina is the top sequencing platform. Types:
MiniSeq (cheap)
MiSeq (Bench top)
MiSeqDX (high throughput)
NextSeq 500 benchtop
HiSeq 2500/3000/4000 - workhorses of sequencing - large scale genomics, production scale.
HiSeq X 5/10
NovaSeq 5000 (for counting)
NovaSeq 6000 (high seuqence coverage)

PacBio sequencers offer longer reads than Illumina. 

**Single read vs. paired end reads:**
* single read = determines the DNA sequence of just one end of each DNA fragment
* paired end = sequence both ends of each DNA fragment making pairing and directionality information available. 
	* Disadvantages: (i) 20% more expensive (ii) measures same fragment twice, thus at same genomic coverage will utilise only half as many unique fragments as single end sequencing. 
* for detecting de novo transcriptome assembly in humans need 100-200 x10^6 paired end reads.
![enter image description here](https://www.yourgenome.org/sites/default/files/images/illustrations/bioinformatics_single-end_pair-end_reads_yourgenome.png)

## Compress & Uncompress Files
A "compressed file" is a single file reduced in size. The name may end with .gz , .bz , .bz2 , or .zip 
A "compressed archive" is a folder containing multiple files combined, and then compressed. The name may end with .tar.gz , .tar.bz2

There are programmes to compress & uncompress files:
- ZIP , extension .zip , program names are zip/unzip
- GZIP extension .gz , program names are gzip/gunzip
- BZIP2 extension .bz/.bz2 , program names are bzip2/bunzip2
- XZ extension .xz . A more recent invention. Technically bzip2 and (more so) xz are improvements over gzip , but gzip is still quite prevalent for historical reasons and because it requires less memory.
- [BGZIP extension](http://www.htslib.org/doc/tabix.html). Specific for bioinformatics. Allows random access to content of compressed file. Can decompress with gzip but only bgzip creates a bgzip file.

To compress a file:
1. efetch sequence file usually in fasta .fa format `efetch sequence id`
2. compress file name with gzip `gzip sequence_name.fa`

You can read compressed files without uncompressing using `zcat` on Linux (`gzcat` on macOS)
`zcat sequence_name.fa.gz | head`

To uncompress file: `gunzip sequence_name.fa.gz` creates the file sequence_name.fa

To compress multiple files together us `tar` Tape Archive:
`tar czfv sequences.tar.gz AF086833.fa AF086833.gb`
Means that we want to create `c` , a compressed `z` , file `f` , in verbose `v` mode in a file called `sequences.tar.gz`. The 2 files to add are listed at the end.

The best was to compress multiple files (if there are many) is to put all files into a folder and then compress that entire directory:
`mkdir sequences`
`mv sequence_names.* sequences/`
`tar czvf sequences.tar.gz sequences/*`

rsyncable archive:
- `rsync` tool synchronises files by sending only the differences between existing files. 
- When files are compressed, you cant use `rsync` but you can use `gzip --rsyncable` flag. This allows gziped files to by synced much faster. e.g. `tar -c sequences/* | gzip --rsyncable > file.tar.gz`

fastq.gz  = compressed version of fast file (needs unzipping before analysing)

[NCBI Format Guide](https://www.ncbi.nlm.nih.gov/books/NBK242622/)

# Experimental Design
 
**Variability** in results:
* need **replicates** to capture breadth of isolate noise
* Technical replicates = repeat library preparations from the same RNA sample —> avoid batch effects & lane effects. Should multiplex same sample over different lanes of same flowcell.
* Biological replicates = parallel measurements on different samples i.e. RNA from independent cells/tissues. Most RNA-seq have 3 biological replicates but ideally need 6 per condition to improve statistical power.
 
**Artificial RNA spike-ins**
* used to accurately quantify absolute transcript concentration. 
* RNA of known quantities is used for calibration eg ERCC. R package `erccdashboard`. Different spike-in controls are needed for each RNA type.
* dont use spike ins to normalise between different samples (they dont account for differences in amount of starting material).
 
**Blocking and randomise**
* Randomly choose which samples to treat and sample
* Block samples into groups based on known sources of variation (sex, weight, cell cycle status) - subexperiments in each block increases sensitivity.
 
# Sequencing Data

Sequencing data are stored in GenBank, FASTA & FASTQ formats. Whilst GenBank & FASTA are usually curated, FASTQ is experimentally obtained.

**1. GenBank format**

- oldest, designed to be human readable --> therefore not optimised for data analysis.
- fixed-width format. first 10 characters form an identifier column.
- convert to simpler format using ReadSeq.
- RefSeq (NCBI Reference Sequence project) database provides reference standards from all the data in GenBank. RefSeq recods have accession numbers that begin with two capitals followed by underscore eg NP_  

**2. FASTA format**

- Header: `>` symbol on the FASTA header line indicates a FASTA record start. The header line may contain an arbitrary amount of text (including spaces) on the same line. This text may include structured information e.g. accession numbers. 
- A string of letters - the sequence ID follows the `>` symbol according to [IUPAC nucleotide guidelines](https://www.bioinformatics.org/sms/iupac.html).
	- `ATGC` represent nucleotides.
	- `N` indicates the base could be any of `ATGC`
	- `W` indicates either of `A` or `T`
	- gaps are represented by `.` or `-`
- sequence lines should be short
- tools can accept letters beyond the IUPAC guidelines and wont be recognised.
- sequence lines should wrap at the same width (exception being the last line).
- use capital letters (lower case letters were previously used to represent repetitive regions) as some tools skip over lower-case regions.

**3. FASTQ format**

- format by which all sequencing instuments represent data.
- each sequence base is associated with a quality score.
-  end in **SRR** and **SRX**
- FASTQ files are uncompressed & large. 

FASTQ format consisted of 4 sections:
1. **Header**: instead of `>` in FASTA, FASTQ start with `@`. Followed by the sequence ID with option for more text which  contains:
	- `@<machine_id>:<run ID>:<flowcell ID>:<flowcell lane>:<tile in the lane>:<x-pos>:<y-pos of the cluster in the tile><read single or paired-end>:<is filtered Y or N>:<control number 0 if none>:<index sequence>`
2. **Sequence**: usually, but not necessarily, on a single line!
3. Starts with `+` followed by optional text  e.g. sequence ID again
4. **Quality scores** of each base in the Sequence (from section 2) - must be same length & wrapped in the same way as section 2. Each character represents a numerical value from the [Phred Score, Q](https://en.wikipedia.org/wiki/Phred_quality_score).

FASTQ example:
5. `@SEQ_ID`
6. `GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT`
7. `+`
8. `!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65 `

![enter image description here](https://www.researchgate.net/profile/Morteza_Hosseini17/publication/309134977/figure/fig2/AS:417452136648711@1476539753452/A-sample-of-the-FASTQ-file.png)

* **Base calling** = deduce the nucleotide letter code sequence from the fluorescence signal edited when incorporated into the sequence read. Imperfect. Information on [base calling here](https://academic.oup.com/bib/article/12/5/489/268399).
* **Phred score, Q** = proportional to probability that a base call is incorrect. P = 10 ^ (-Q/10)
	*  10 = 1 in 10 bases are wrong (90% accuracy); 
	* 20 = 1 in 100 bases are wrong (99% accuracy). 
	* Higher Phred = higher quality. 
	* Maximum score is 40.  
	* Each Phred score is represented as a single ASCII character but may represent a 2 digit score, e.g. score I = Phred of 40. 
* Using the current standard (Sanger +33 format), ASCII character symbols are mapped to scores as follows:

ASCII character = !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
Phred score ----= 0------5---10---15-----20-----25---30------35-----40
Quality----------= worst-----------------------------------------------best

* There are different versions of FASTQ quality score encoding. Depends on the sequence technology and the base caller assignment used (eg Bustard, RTA0 HiSeq X). The newer Sanger +64 format has the following mapping:

ASCII character = @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi
Phred score ----= 0-------5-----10-------15-----20-------25---30----35----40
Quality----------= worst-------------------------------------------------------best

- Note that in Sanger +64, characters `@` `A` `B` `C` etc map to low quality (score <5) whereas in Sanger +33 these same characters map to high quality (>30).

As a general rule remember that:
- !"#$%&'()*+,-.  = low quality 1/10
- -./0123456789 = medium quality 1/100
- ABCDEFGHI = high quality 1/1000

Easily distinguish between Sanger +33 and +64 by the presence of:
- lower case letters = +33 version. 
- 0 character = +33 version.  When the quality scores do not contain 0, it is either Solexa +64, Illumina 1.3+ Phred+64, Illumina 1.5+ Phred+64.
- characters J - Z = +64 version

Depending on the tool being used you can either (a) account for different version or (b) convert data to correct coding, e.g. using `seqtk` 
- e.g. `seqtk seq -Q64 input.fq >output.fa`
- converting Illumina FASTQ file 1.3 (Phred + 64) > version 1.8 (Phred +33) can also use: 
`sed -e '4~4y/ @ABCDEFGHIJKLMNOPQRSTUVWXYZ [\\]^_abcdefghi/!"#$%& '\ ' '()*+ , -.\/0123456789:; <= >? @ABCDEFGHIJ /' originalFile.fastq`

## Sequencing Database Repositories:
**[Sequencing Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)**: main repository for raw sequencing reads from next generation sequencing. Includes NCBI & European Bioinformatics Institute & DNA Databank of Japan
**[European Nucleotide Archive, (ENA)](http://www.ebi.ac.uk/ena)**: similar to SRA except files are stored directly as FASTQ & BAM formats (dont need to convert with `sratoolkit`)
**[Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/)**: All GEO sequences are in SRA (both NCBI). GEO consists of functional genomic data repository including RNA-Seq, ChIP-seq, RIP-seq, HiC-seq, methyl-seq etc.  Only stores gene expression level results - the sequence data is deposited in SRA. (i.e. projects can have two locations for data). GEO accession ID start with GSE
**GenBank**: everything else. Divisions include Genomes (complete genome assemblies) & WGS (whole genome shotgun)
**[Transcriptome Shotgun Assembly (TSA)](https://www.ncbi.nlm.nih.gov/genbank/tsa/)**: Transcriptome assemblies. Division of GenBank 
**[ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)**: in addition to GEO is the main repository for functional genomics data e.g. RNA-seq. From EMBL-EBI.

## Accessing Sequence Data

**NCBI & Entrez**
- ginormous data store of sequence data.
- Entrez is NCBIs primary text search integrating PubMed with 39 other databases
- Entrez enables you to access NCBI databases via the command line. 
- Accession number applies to the complete database record for that entity. Updates and revisions remain under the same accession number but with a different version number.
- Prior to 2015 NCBI versions began with GI.
- Entrez web API allows us to query NCBI data using the URL construct: 

`https://service.nih.gov/?param1=value1&param2=value2&param3=value3`
- place `\` before each `&` to escape the command line `&` meaning or alternatively place ' ' around the whole URL.
- For example to access AF086833.2 in FASTA format via the command line: 

`curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=AF086833.2&db=nuccore&rettype=fasta'`

- Entrez direct tool simplifies the command line access e.g. `efetch`, `esearch`, `xtract`

`efetch -db=nuccore -format=gb -id=AF086833 | head`

#e.g. accession number AF086833 in Genbank format. `efetch -db=nuccore -format=gb -id=AF086833 > AF086833.gb`
#e.g. accession number AF086833 in Fasta format. `efetch -db=nuccore -format=fasta -id=AF086833 > AF086833.fa`
#efetch can take additional parameters and select a section of the sequence. `efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=3`
#It can even produce the sequence from reverse strands: `efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=1`
>gb|AF086833.2|:1-5 Ebola virus - Mayinga, Zaire, 1976, complete genome CGGAC

`efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=2`

>gb|AF086833.2|:c5-1 Ebola virus - Mayinga, Zaire, 1976, complete genome GTCCG

Note strand 2 is the reverse complement (not simple complement) represented by `c5-1` tag.

- searching Entrez Direct:
	- use the project accession number from the published paper to search the data that comes with it

 `esearch -help`
 `esearch -db nucleotide -query accession_number`

To fetch the data from a search pipe it to efetch: `esearch -db nucleotide -query accession_number | efectch -format=fasta >genome.fa` This generates a XML file.

To navigate and select parts of the XML file from `efetch` pipe it to `xtract`
`efetch -db taxonomy -id 9606,7227,10090 -format xml | xtract -Pattern Taxon -first TaxId ScientificName GenbankCommonName Division`

**[ENSEMBL](http://www.ensembl.org/index.html)** : has the most detailed annotations of genomes. Download the latest Ensembl release of the human transcriptome `curl -0 ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz`

**UCSC Table Browser**: Stores large scale result tiles produced by consortia e.g. ENCODE. Contains multiple sequence alignment datasets across entire genomes.
* UCSC and Ensembl use different naming conventions (which impacts on analyses) - try to stick to one.
* must use Ensembl FASTA file with Ensemble GTF file. You cannot mix Ensembl with UCSC without modifying first and this isn't advised as they have different scaffolds.

## Downloading from Sequence Database Repositories

`ml ncbi-vdb`
`ml fastq-tools`
`ml seqtk`

1. Ctrl & F in manuscript PDF to find unique ID accession numbers. Summarised [here](https://www.ncbi.nlm.nih.gov/guide/howto/submit-sequence-data/).
2. For `GSE****` accession numbers go to: https://www.ncbi.nlm.nih.gov/geo/ and search the `GSE****` ID
3. On the [GEO Accession Display](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98290) page click the BioProject or BioSample ID under Relations `PRJNA****`
	- NCBI BioProjects accession IDs start with `PRJNA****`
	- NCBI BioSample start with `SAMN****` and `SRS****`. They describe biological source material
4. On the [BioProject page](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA384592) click on SRA Experiments links under Project Data
   -  `esearch`  the BioProject number: `esearch -db sra -query PRJNA******` . The reveals how many sequencing runs there are for the project ID.
   - Format the output as a RunInfo.csv format: `esearch -db sra -query PRJNA****** | efetch -format runinfo > info.csv`
   - To see SRR IDs, head the first column: `cat info.csv | cut -f 1 -d ',' | head`
	- Filter for `SRR` & store the first 10 IDs: `cat info.csv | cut -f 1 -d ',' | grep SRR | head > ids.txt`
6. On the [SRA page](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=384592) show 200 results & then click on Send Results to Run Selector. This creates the SRA Run Selector table which has useful information on each sequencing file.
	- [SRA Handbook](https://www.ncbi.nlm.nih.gov/books/NBK47528/)
	- SRA experiment IDs start with `SRS****`. Unique sequencing library for a sample.
	- SRA run IDs start with  `SRR****` and `ERR****`. Data file linked to the sequencing library.
7. On the [SRA Run Selector page](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_2827488_130.14.18.97_5555_1539957447_3051893689_0MetA0_S_HStore&query_key=3) click Accession List icon under Download. This downloads a text file of each Sequencing run ID.
8. In terminal download each of the `SRR****` IDs in the txt file using the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) tool `fastq-dump` - this converts the SRA sequences to a FASTQ file. 
	- in command line type `fastq-dump SRR1553607`. Generates `SRR1553607.fastq` file

The `while loop` allows you to do this for all SRR IDs in a single command and `sbatch` speeds up the process by parallelising the request: 
```bash
while read line; 
do 
	sbatch -N 1 -c 1 --mem 32 --wrap="fastq-dump --split-files --accession $line"; 
done < SRR_Acc_List.txt
```
`N` is the number of nodes (usually leave at 1)
`c` is the number of cores (i.e. threads - so you could change the STAR command to `runThreadN 8` with this example)
`--mem` is the amount of memory you want
Can remove the `slurm` files after download (that is just the log): `rm slurm-*`

Alternatively using the `ids.txt` file list of SRR run IDs & invoke `fastq-dump` on each ID: `fastq-fump -X 10000 --split-files SRR******` but this can be done in one go using: `cat ids.txt | xargs -n 1 fastq-dump -X 10000 --split-files $1`

- **Paired-end reads** are concatenated by SRA (because that is how the instrument measures them). Therefore need to separate these into 2 different FASTQ files: 1 forward read; 1 backward read using the command: `fastq-dump --split-files SRR1553607`. Generates files `SRR1553607_1.fastq` & `SRR1553607_2.fastq`
	* know the origin of each read (forward vs reverse) - encoded in read name - some analysis tools require combining the 2 files into 1. The **forward** read is **filename_1** and **backward** read is **filename_2**
	* Need to process the read as Split Reads/files

**Alternative approach to downloads:**
- Copy the URL Address for the relevant FASTQ files from the column in the [SRA run selector table](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_110487706_130.14.18.97_5555_1538386558_3578612985_0MetA0_S_HStore&query_key=2). 
- in terminal change to the target directory `cd`, then in command line use `wget` or `curl` tools to download e.g. `curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972739/SRR
1972739.sra` or  	`wget link_copied_from the_website`
- These generate SRA files which need to be converted to FASTQ using: `fastq-dump -X 10000 --split-files SRR1972739.sra` # we specified to download the first 10,000 reads using `-X 10000` and split the "spots" using `--split-files`. Spots are reads.
- If there are many samples then download summary (right click on TEXT) & copy link location: then in command line run: `wget -O samples_at_ENA .txt "<paste LINK >"` #the quotation marks are crucial 
- use 11th column of TEXT file (Fastq file top) to feed the URLs of different samples `cut -f11 samples_at_ENA . txt | xargs wget`

## View & Interrogate the downloaded data: 
 - `more file_name.README` Press `space` to move forward; `b` to move back, `q` or `ESC` to exit.
 - print the content of the file to the terminal window: `cat file_name`
	 - scan the text to find columns and rows of interest
- count lines, words & characters: `cat file_name | wc`
- how does the file start: `cat file_name | head`
- find information on a specific gene or transcript by searching for text matching the gene name: ` cat file_name | grep gene_x`
	- `grep -v` will show lines that don't match
	- `cat file_name | cut -f 2 | grep gene_feature | head` to select genes, this command shows the gene_feature in column 2
	- `cat file_name | cut -f 2 | grep gene_feature | wc -l` counts the number of genes present
- if you plan to use this grepped data to interrogate further then place it in a separate file to analyse: `cat file_name | cut -f 2 > new_file.txt`
	- sort into identical consecutive entries: `cat new_file.text | sort`
	- collapse duplicated identical words: `cat new_file.txt | sort | uniq | head`
	- print counts: `cat new_file.txt | sort | uniq -c | head`
	- unique types of features: ` cat file_name | cut -f 2 | sort | uniq -c | wc -l`

The [`sra-stat` program](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-stat) can rapidly generate an XML report on the sequence data: 
`sra-stat --xml --quick SRR1553610`	
`sra-stat --xml --statistics SRR4237168` generates read length statistics.

## Manipulating FASTQ & FASTA files (Optional Step)
`ml seqtk` 
Installed [SeqKit](https://github.com/shenwei356/seqkit) & [csvtk](https://github.com/shenwei356/csvtk)

Overview of FASTQ data: `seqkit stat *.gz` or `seqkit stat *.fastq`
 
Convert FASTQ file into 3-column tabular format (1st col name/ID, 2nd col sequence, 3rd col quality) then add optional columns e.g. sequence length, GC content: `seqkit fx2tab --name --only-id --gc *.gz`

Extract a subset of sequences from a FASTQ file: ` seqkit sample --proportion 0.001 *.gz | seqkit seq --name --only-id > id.txt`
ID list file: `head id.txt` prints the relevant run IDs with specified proportion 0.1%
Search by ID list file: `seqkit grep --pattern-file id.txt *.gz > subset.fq.gz`
 
 Find degenerate bases in FASTQ file:  `seqkit fx2tab` converts FASTQ to tabular format & outputs the sequence in a new column.  `seqkit fx2tab --name --only-id --alphabet *.gz | csvtk --no-header-row --tabs grep --fields 4 --use-regexp --ignore-case --pattern "[^ACGT]" | c svtk -H -t cut -f 1 > id2.txt`
 
 Use grep text search to filter the table: `seqkit grep --pattern-file id2.txt --invert-match *.gz > clean.fa`
Or locate degenerate bases K & N: `seqkit grep --pattern-file id2.txt viral.1.1.genomic.fna.gz | seqkit locate --ignore-case --only-positive-strand --pattern K+ --pattern N+`

Remove duplicated sequences: `seqkit rmdup --by-seq --ignore-case *.fq.gz > uniq.fq.gz`

Locate a specific motif in FASTQ format: `seqkit locate --degenerate --ignore-case --pattern-file *.gz` use the degenerate flag to identify pattern of interest
 
## Quality Control (QC)
Quality control (abbreviated as QC from now on) is the process of improving data by removing identifiable errors from it. QC is performed at different stages:
- Pre-alignment: “raw data” - the protocols are the same regardless of what analysis will follow
	- **FastQC** on raw sequenced reads
	- reliability of sequencing decreases along the read. Trimming works from the end of the read & removes low quality measurements.
	- trim adapters
- Post-alignment: “data filtering” - the protocols are specific to the analysis that is being performed.

### FastQC
FastQC on raw reads FASTQ file using the [FastQC program](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - need to install. FastQC generates its reports by evaluating a small subset of the data and extrapolating those findings to the entirety of the dataset. Many of the metrics are only computed on the first 200,000 measurements then are being tracked through the rest of the data.

Check for:
* PCR duplicates
* adapter contamination
* rRNA + tRNA reads
* unmappable reads (contaminating nucleic acids) - FASTQC doesn’t check this
* each test will either = pass; warn; or fail. Fail is expected in some cases and does not mean the sequencing needs repeated.
 
__FastQC Workflow__
`ml FastQC`
`ml python`
`ml pip`

make a folder to store the results: `mkdir fastqc_results`
run FastQC on each sequencing file: `fastqc SRR5* -o fastqc_results`
$ for SAMPLE in WT_1 WT_2 WT_3 WT_25 # random selection of samples: `mkdir fastqc_results/${SAMPLE}`
if there are many sequencing files can speed up by running as `sbatch`

Look at the results: `ls fastqc_results/ERR458493_fastqc/`
Summarise multiple FastQC outputs using the [MultiQC tool](http://multiqc.info/):
[run `multiqc`](http://multiqc.info/docs/#running-multiqc) within the `fastqc_results` folder and use the folder names (WT_1 etc.) as prefixes to the sample names in the final output: 
Go to the folder with the fastqc files in and simply run: `multiqc .`
Open the MultiQC html report: `open multiqc_report.html` or open in the browser
[Interpret](https://www.youtube.com/watch?v=qPbIlO_KWN0) the MultiQC report: general stats, quality scores (Phred score), 

### Trim low quality bases & adapters & rRNA
There are many QC tools available (most in bash but some in R - bioconductor) each with basic QC methods plus unique functionality. Best ones include:
`Trimmomatic`  application note in Nucleic Acid Research, 2012, web server issue
`BBDuk` part of the BBMap package
`FlexBar` Flexible barcode and adapter removal  published in Biology, 2012
`CutAdapt`  application note in Embnet Journal, 2011 - advised by Nobby
`Trimer Galore` is a wrapper around `Cutadapt` - this performs quality trimming then adaptor trimming in one step.

`ml cutadapt`
`ml trimomatic`
`ml Trim_Galore`
`ml FASTX-Toolkit`
`ml Bowtie2`

`trim_galore -q 20 --gzip -o trim_galore SRR5*_1.fastq`

for multiple sequences can parallelise by using for loop & sbatch:

```bash
for file in ~/working/oliver/projects/airals/fastq_files/D7_samples/SRR5*_1.fastq
do
	sbatch -N 1 -c 1 --mem 8 --wrap="trim_galore -q 20 --length 20 --gzip -o trimmed_D7 $file";
done
```
@Raphaelle used:
adapter removal: `fastx_clipper -Q 33 -l 24 -a $adapt -i ${paths}${file} > ${paths}${data}-clipped.fastq`

rRNA removal using bowtie (the indexed fasta sequence of ribosomal RNA as obtained from UCSC): `bowtie -q -p 8 -v 0 --un ${paths}${data}-no_rRNA.fq $ribosomal ${paths}${data}-clipped.fastq > ${paths}${data}-rRNA.sam`
 
### Error Correction (Optional step)
Nobby doesn't routinely do this. 
Sequencing instruments occasionally make random errors. Can recognise errors by computing k-mer density. K-mers are all possible short subsequences of length `k` contained in a string (sequence) that occur many times. They can be 2-5 bases long. K-mers that contain errors will be rare so can be used in error correction but also genome identification/classification & alignement.

- The 2 base long k-mers ( 2-mers ) are AT , TG , GC and CA
- The 3 base long k-mers ( 3-mers ) are ATG , TGC and GCA
- The 4 base long k-mers ( 4-mers ) are ATGC , TGCA
- The 5 base long k-mer ( 5-mer ) is ATGCA

FaASTQ error correction programs correct or remove reads that appear to have errors in. Very useful when data is of high coverage.

`ml BBMap`
`bbmap` package using `tadpole.sh` error corrector
`tadpole.sh in=SRR5*_1.fastq out=tadpole.fq mode=correct`

### Sequence patterns
Sequence patterns can be summarised probabilistically into motifs e.g. where `GC` is followed by `A` 80% of time.
Adapters are very simple sequence patterns that can be identified by a pattern. 
When sequences are on a single line (line orinentated), `grep` with `--colour=always` allows us to identify the mataches e.g. `cat SRR5483788_1.fastq | grep --color=always CTGTCTCTTATACA | head -2`

Searching sequences for patterns:
`grep` and extended grep `egrep` work on single line.
`dreg` and `fuzznuc` Emboss tools wra over many lines

Regular expressions:
`^` matches the beginning of the line
`.` matches any character
`{m,n}` matches the preceding elements at least m but not more than n times.

## Sequence Alignment
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

### Types of alignment

- **Global alignment** is where every base of both sequences has to align to another matching base, to another mismatching base, or to a gap in the other sequence.
- **Local alignment** algorithms look for the highest scoring subregions (or single region). Local alignments are used when we need to find the region of maximal similarity between two sequences.
- **Semi-global alignment** (global-local, global) is a method mixture between the global and local alignments. It attempts to fully align a shorter sequence against a longer one. They are used when matching sequencing reads produced by sequencing instruments against reference genomes. Majority of data analysis protocols rely on semi-global alignment.
- Multiple sequence alignments uses 3 or more sequences (i.e. not pairwise alignment). 

### Short read aligners
In 2005 high throughput short-read sequencers changed the face of sequencing which became much cheaper and accessible. Previously sequencing was laborious and expensive and the focus was on producing accurate long reads (1000bp). Sequencing reads longer improves alignment.  Sequencers now produce millions of short reads (50-300bp). Aligners have thus changed to adapt to short reads - rapidly select the best mapping at the expense of not investigating all potential alternative alignments.
- Short read aligners are incredible! They match > 10,000 sequences per second against 3 billion bases of the human genome.
- There is large variation in results between different aligners. A tool that prioritises finding exact matches will do so at the exense of missing locations and overlaps and vice versa.
- Limitations:
	- finds alignments that are reasonably similar but not exact (algorithm will not search beyond a defined matching threshold)
	- cannot handle very long reads or very short reads (<30bp) (become inefficient)

### Choosing an Aligner (STAR, bwa, bowtie2, TopHat2)
When choosing an aligner, we need to decide what features the aligner should provide:
-   Alignment algorithm: global, local, or semi-global?
-   Is there a need to report non-linear arrangements?
-   How will the aligner handle INDELs (insertions/deletions)?
-   Can the aligner skip (or splice) over large regions?
-   Can the aligner filter alignments to suit our needs?
-   Will the aligner find chimeric alignments?

For RNA-seq we need to:
- align reads to the reference transcriptome index (required transcripts to be known and annotated in the reference)
- identify novel splice events (using reads that cant be aligned to reference transcriptome)

![enter image description here](https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/resize/user/18/Figure19-700x527.png)

* Multiple alignment programmes available for RNA-seq, each specialising in detecting different factors eg structural variants; fusion transcripts
* Straight forward RNA seq for differential gene expression analysis = use STAR [STAR manual PDF](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
	* Efficient
	* Sensitive
	* But large number of novel splice sites (caution)
* **TopHat** = popular aligner – wrapper around the genomic aligner Bowtie

### Reference Genomes
Reference genome sequence repositories:
**GenCODE** - Generally for human data use 
[ENSEMBL](http://www.ensembl.org/index.html) - use for most other species
[UCSC](https://genome.ucsc.edu/) - comprehensive comparative genomics data across genomes         	https://genome.ucsc.edu/cgi-bin/hgTables
ENCODE
iGenomes
NCBI
Mouse Genome Project
Berkeley Drosphilia Project
RefSeq
[RNA-Central](http://rnacentral.org/) - non coding RNA sequencing database

__Reference Genome File Formats__
Reference sequences are **FASTA files.** Reference sequences are long strings of ATCGN letters.  File formats store start sites, exon, introns. One line per genomic feature.

Reference annotations are **GFF** or **GTF** format.

GFF = General Feature Format
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
 
**GTF** = Gene Transfer Format. More strict than GFF. Same as GFF2 but 9th field expanded into attributes (like GFF3). http://mblab.wustl.edu/GTF2.html

**BED Format** is the simplest annotation store
3 compulsory fields: chromosome & start & end.
9 optional fields: name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        	field number can vary from 3 - 12. Must be consistent within a file.
indicates region with 0-based start and 1-based end position (GFF & GTF are 1-based in both start and end) Aligning Reads

### Download Reference Sequence (FASTA) & Annotation (GTF) files
Keep reference genomes in a more general location rather than the local folder, since you will reuse them frequently: `/home/camp/ziffo/working/oliver/genomes/sequences` and then make a new directory `mkdir`

**GENCODE process:**
https://www.gencodegenes.org/releases/current.html
find the latest human reference genome: currently this is Release 28 (GRCh38.p12)
Copy the Link address to download the FASTA file to the GENOME SEQUENCE  PRI (primary) - bottom of 2nd table. 
Then in command line download to the appropriate directory: `wget  [paste link address]` then `gunzip filename`
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

Send cmd to generate index as `batch job` to cluster:
`sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runMode genomeGenerate --genomeDir /home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index --genomeFastaFiles /home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf --sjdbOverhang 59 --runThreadN 8"`

Check it is running using `myq`
If it isn’t seen there then `ls` the folder → look for output file named “slurm…”
Open slurm file: `more slurm…` which will explain outcome of file eg FATAL INPUT PARAMETER ERROR
 
To make the command shorter and easier to interpret you can assign names to files using $name. E.g
`REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa`
`ANO=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28.primary_assembly.annotation.gtf`

`sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runMode genomeGenerate --genomeDir /home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index --genomeFastaFiles $REF --sjdbGTFfile $ANO --sjdbOverhang 59 --runThreadN 8"`

Alternative approach is to build Index using [snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) which parallelises the indexing process:

 - Write the Snakefile in Atom and save it in the relevant directory near the index. Define input, output and shell script.

`input: ref="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa", starref="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR",
output: "/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index"+"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR"+"/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index/chrName.txt"
shell: 
"STAR --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.ref} --sjdbOverhang 59 --runThreadN 8"`

 
### 2. Align each FASTQ file to the Index
ml Python/3.5.2-foss-2016b
ml SAMtools

* Sample distributed over X flow cell lanes → X fastq files per sample
* STAR merges the X files if multiple file names are indicated (other align tools don't)
* Separate the file names with a comma (no spaces)
* Create directory to store STAR output `mkdir alignment_STAR`

By allocating all file names to the `$INPUT` term it means all the FASTQ files are read together (see example further down) but it is best to do each separately and use the [snakemake command](http://slides.com/johanneskoester/snakemake-tutorial#/) - this makes it easier to see if there is an alignment error in each individual sequencing file:

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
```

dry run workflow:
`snakemake -np /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam`
execute workflow:
`snakemake /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam`

`snakemake /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam`

* List fast.qz files separated by commas and remove white spaces:

`INPUT= ls -m /home/camp/ziffo/working/oliver/projects/airals/fastq_files/SRR5*.fastq| sed 's/ //g' | echo $INPUT | sed 's/ //g'`

`ls -m` list, fill width with a comma separated list of entries all the fastq files.
`sed` remove a space from each file name
`echo` displays a line of text - in this case it lists all the file names  - in this case the pipe to echo is to remove new spaces that are created between multiple lines.
no space after FILES - with space after it thinks FILES is command.
`sed` = stream editor - modify each line of a file by replacing specified parts of the line. Makes basic text changes to a file - `s/input/output/g`

`
sbatch -N 1 -c 8 --mem 40G --wrap="STAR --runThreadN 1 --genomeDir /home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index --readFilesIn /home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/SRR5483788_1.fastq,SRR5483789_1.fastq,SRR5483790_1.fastq,SRR5483794_1.fastq,SRR5483796_1.fastq,SRR5483795_1.fastq --outFileNamePrefix /home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/D7_samples/D7_ --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic"`

Optional additions in STAR manual: e.g.
if files are compressed add `--readFilesCommand gunzip -c`

Check the summary of the output:
`cat FILENAME_Log.final.out`
![summary of alignment](https://user-images.githubusercontent.com/33317454/37438378-095cb442-282d-11e8-95fd-bfecefae5b75.png)

The [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  has all the explanations on how to write the STAR command and fine tune parameters including:
* multi-mapped reads
* optimise for small genomes
* define min & max intron sizes 
* handle genomes with >5000 scaffolds
* detect chimeric & circular transcripts

**STAR output files:**
 - Aligned.sortedByCoord.out.bam - the loci of each read & sequence
 - Log.final.out - summary of alignment statistics
 - Log.out - commands, parameters, files used
 - Log.progress.out - elapsed time
 - SJ.out.tab - loci where splice junctions were detected & read number overlapping them
 - Unmapped.out.mate1 - fastq file with unmapped reads

More information on these is in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) chapter 4 page 10.

**ENCODE options:**
`outFilterMultimapNmax 1` max number of multiple alignments allowed for a read - if exceeded then read is considered unmapped i.e. when 1 is used this only identifies unique reads and removes multimapped reads. This is generally accepted.

#`STAR` will perform the alignment, then extract novel junctions which will be inserted into the genome index which will then be used to re-align all reads
#`runThreadN` can be increased if sufficient computational power is available

### 4. [Indexing read alignments](http://software.broadinstitute.org/software/igv/bam)
Generate the SAM file > convert to BAM format > index the BAM file.
Use `samtools` package to index. Load `samtools` in command line: `ml SAMtools`
Create BAM.BAI file with every BAM file to quickly access the BAM files without having to load them to memory.
Run samtools index cmd for each BAM file once mapping is complete:
`samtools index alignment_STARAligned.sortedByCoord.out.bam`

`samtools sort alignment_STARAligned.sortedByCoord.out.bam > STARAligned.sortedByCoord.out.sam`

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

Look at the Directed Acyclic Graph (DAG):
`snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg`

@Raphaelle used:
Map with tophat2 (however please note that I did this step in 2016 and if I would need to restart now, I would rather use STAR) : `tophat -g  1 -p 8 -G $geneModel --library-type fr-firststrand -o ${out}${data} $genome ${paths}${data}-no_rRNA.fq`

## Sequence Alignment Maps (SAM)

SAM files are generic nucleotide alignment format describing alignment of sequenced reads to a reference.. SAM format are TAB-delimited line-orientated (each row represents a single read alignment) text consisting of:
1. Header: meta-data
2. Alignment: longer with information on alignment.

-   The  [SAM format specification](http://samtools.github.io/hts-specs/SAMv1.pdf)  is the official specification of the SAM format.
-   The  [SAM tag specification](http://samtools.github.io/hts-specs/SAMtags.pdf)  is the official specification of the SAM tags.

* Each SAM file has 11 mandatory columns. Despite this aligners vary in how much information they put into these columns - impossible to transform one aligners SAM file into another aligners SAM file.

![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/bam_structure.png)

**SAM header section**
- Header rows start with `@`then tag: value pairs 2 letter abbreviations e.g. `SQ` = sequence directory listing  `SN` sequence name aligned against (from FASTA), `LN`  sequence length, `PG` program version that was ran.
* 1 line per chromosome
* To retrieve only the SAM header `samtools view -H`
* To retrieve both the header & alignment sections `samtools view -h`
* The default `samtools view` will not show the header section

`samtools view -H WT_1_Aligned.sortedByCoord.out.bam`

**SAM alignment section**
* Each line corresponds to one sequenced read
* 11 mandatory columns in specified order:

`<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>`

* If info is unavailable then value can be 0 or * 
* Optional fields follow the mandatory fields

![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/sam_fields.png)

*2nd column FLAG field:*
* stores info on the respective read alignment in one single decimal number
* decimal is the sum of all the answers to Yes/No questions:
![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/sam_flag.png)
* To convert the FLAG integer into plain english [click here](https://broadinstitute.github.io/picard/explain-flags.html).
* When a read has multiple alignments matching equally well the aligner designates one as the primary alignment.

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

__Read Group__
Key feature of SAM/BAM format is ability to label individual reads with readgroup tags. This allows pooling results of multiple experiments into a single BAM dataset to simplify downstream logistics into 1 dataset.
Downstream analysis tools eg `variant callers` recognise readgroup data and output results.

[BAM readgroups GATK website](https://gatkforums.broadinstitute.org/gatk/discussion/1317/collected-faqs-about-bam-files)

Most important readgroup tags: `ID`, `SM`, `LB`, `PL`
![enter image description here](https://galaxyproject.github.io/training-material/topics/introduction/images/rg.png)

### BAM files
- BAM files are binary, compressed and sorted representation of SAM information. 
- They are usually sorted by alignment coordinate which allows fast query on the information by location.
- Can rarely be sorted by read name when read pairs (with identical names) need to be accessed quickly (stored in adjacent lines)
- Used to exchange data as it quicker than SAM.

### CRAM format
Similar to BAM (binary compressed) but smaller as some compression is in the reference genome.
Sometimes you need the reference genome information so these arnt always appropriate.
Supported by samtools
Concerted effort to move from BAM to CRAM.

### 5. Store Aligned Reads as SAM/BAM files
SAM & BAM files contain the same information but in different formats

 - `sam-dump` tool downloads SRA data in SAM alignment format.

**Converting BAM files to and from SAM files**
BAM --> SAM:
`samtools view -h WT_1_Aligned.sortedByCoord.out.bam > WT_1_Aligned.sortedByCoord.out.sam`
Convert a BAM file into a human readable SAM file (including the header): `samtools view -h FILENAME.bam > FILENAME.sam`

Compress a SAM file into BAM format (-Sb = -S -b)" `samtools view -Sb FILENAME.sam > FILENAME.bam`

To peak into a SAM or BAM file: `samtools view FILENAME.bam | head`

Generate an index for a BAM file (needed for downstream tools): `samtools index FILENAME.bam`


## Manipulating SAM/BAM files
There are 4 major toolsets for processing SAM/BAM files:

- [SAMTools](http://www.htslib.org/) - interact with high throughput sequencing data, manipulate alignments in SAM/BAM format, sort, merge, index, align in per-position format
- [Picard](https://broadinstitute.github.io/picard/) - Java tools for manipulating high throughput sequencing data
- [DeepTools](https://deeptools.readthedocs.io/en/develop/) - visualise, quality control, normalise data from deep-sequencing DNA seq
- [BAMtools](https://github.com/pezmaster31/bamtools/wiki/Tutorial_Toolkit_BamTools-1.0.pdf) - read, write & manipulate BAM genome alingment files

__Processing:__
1. **Filter** using BAM Tools
	-  mapping quality: remove poor alignments - eg remove all alignments below Phred scale quality of 20
	- keep only those which are "properly paired" ie forward is looking at reverse (for paired reads)
	- reference chromosome: e.g. only keep mitochondrial genome alingments
2. **Remove duplicates** with Picard
3. **Clean up** with CleanSam Picard tool
	- fixes alignments that hang off ends of ref sequence
	- sets MAPQ to 0 if read is unmapped

SAMTools help page = `samtools --help`
Usage:   `samtools <command> [options]`

5 key SAMTool commands:
1. Indexing
2. Editing
3. File operations (aligning, converting, merging)
4. Statistics
5. Viewing

For each command there are multiple options `samtools COMMAND -X` 
Create BAM with only **unmapped reads**: `samtools view -h -b -f4 FILENAME.bam > unmapped_reads.bam` 
Create BAM with only **mapped reads**`samtools view -hb -F 4 FILENAME.bam > mapped_reads.bam` 
Create BAM with **mapping quality >= 20**`samtools view -h -b -q 20 FILENAME.bam > high_mapq_reads.bam` 
Create BAM with **uniquely aligned reads** (STAR gives uniquely aligned reads a mapping quality of 255 so you can use samtools to pull all reads with mapping quality = 255 only (using samtools command, option -q, = 255) `samtools view -h -q 255 FILENAME.bam > uniquely_aligned_reads.bam`

- `-h` is used to print the header (always needed).
- sam files are human readable, bam are compressed. BAM are much smaller
	- `-b`will produce a BAM file. `-s` will produce a SAM file.
	- for SAM files you can run other commands on them eg head FILENAME.sam whereas BAM files need to be run through samtools i.e. `samtools view FILENAME.bam | cut -f 2 | head`

Create BAM file with only **reads aligned to reverse strand**:
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

## Quality control of Aligned Reads
Analyses now switch from command line to R studio.

After aligning and before performing downstream analyses check for:
1. Excessive amounts of reads not aligned
2. Obvious biases in the read distributions
3. Similarity between replicate samples
 
 The main STAR output file for downstream analyses are
 **Aligned.sortedByCoord.out.bam** & **Log.final.out**.

**out.mate1 files** are fastx files that are uncompressed and can be made smaller using gzip

__Alignment Assessments__
 
 **1. Check that alignment rate of RNA-seq reads is > 70%**
 
 Check the aligner's output: `cat FILENAME_Log.final.out`
 - most important number = **uniquely mapped reads**
 - if using >2 BAM files then visualise alignment rate for each e.g. using MultiQC in **R studio**: see https://github.com/friedue/course_RNA-seq2015/blob/master/01_Alignment_visualizeSTARresults.pdf
Mount files onto laptop: Right click on Finder --> Connect to server --> Connect to Luscombe Lab
- `infiles <- list.files(path="/Volumes/lab-luscomben/working/oliver/projects/rna_seq_worksheet/alignment_STAR", pattern = "WT_1_Log.final.out", full.names = TRUE)`
- `align.results <- lapply(infiles, function(x) read.table(x, sep="|", strip.white = TRUE, stringsAsFactors = FALSE, skip = 3, fill = TRUE, header = FALSE))`
- `typeof(align.results)`
- `head(align.results[[1]])`
-  `align.results <- lapply(align.results, function(x) transform(x,V2 = as.numeric(gsub("%", "", x$V2) )))`

 **2. Calculate number of alignments in each BAM file**
- easiest to do using a line count ` samtools view Aligned.sortedByCoord.out.bam | wc -l`
- unmapped reads in the BAM file & also multiple instances of the same read mapped to different locations will also be counted (latter only if multi-mapped reads were kept) so run specific tools to indicate FLAG values too.
- `samtools flagstat` assesses the FLAG field and prints a summary report: `samtools flagstat Aligned.sortedByCoord.bam`


### Visualising Data in R

`ml RSeQC`

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

__Visualise the output of  `RSeQC` quality controls__
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

## Bias Identification

**Typical Biases of RNA-seq:**
- many reads aligned to introns indicates: 
	- incomplete poly(A) enrichment 
	- abundant presence of immature transcripts
- many reads align outside of annoted gene sequences (intergenic reads) indicates:
	- genomic DNA contamination
	- abundant non-coding transcripts
- over representation of 3' portions of transcripts indicates RNA degradation

__Read distribution__
- mRNA reads should mostly overlap with exons. Test this with `read_distribution.py` script
	- counts number of reads overlapping with various genes & transcript associated genomic regions (introns and exons)
- download BED file from [UCSC genome](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=685446505_FqnRnlREChczp8SYDIJOSvLwBshv&clade=other&org=S.+cerevisiae&db=sacCer3&hgta_group=genes&hgta_track=sgdGene&hgta_table=0&hgta_regionType=genome&position=chrIV%3A765966-775965&hgta_outputType=primaryTable&hgta_outFileName=).

`read_distribution.py -r sacCer3.bed -i /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/alignment_STAR/WT_1_Aligned.sortedByCoord.out.bam`

Output: 
Total Reads                   937851
Total Tags                    947664
Total Assigned Tags           0

**Group               		Total_bases         Tag_count           Tags/Kb**
CDS_Exons           	8832031             0                   0.00
5'UTR_Exons         0                   0                   0.00
3'UTR_Exons         0                   0                   0.00
Introns             69259               0                   0.00
TSS_up_1kb          2421198             0                   0.00
TSS_up_5kb          3225862             0                   0.00
TSS_up_10kb         3377251             0                   0.00
TES_down_1kb        2073978             0                   0.00
TES_down_5kb        3185496             0                   0.00
TES_down_10kb       3386705             0                   0.00

Visualise this output using this [R script](https://github.com/friedue/course_RNA-seq2015/blob/master/02_Alignment_QC_visualizeReadDistributionsAsBarChart.R).

__Gene body coverage__
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

__mRIN calculation using tin.py__
- RNA integrity number (RIN) is rarely reported in public data repositories.
- Instead determine a measure of mRNA degradation in silico using RSeQCs tin.py script to produce a TIN.
- TIN 0 (worst) - 100 (best). TIN 60 = 60% of transcript has been covered.
- tin.py uses the deviation from an expected uniform read distribution across the gene body as a proxy

`tin.py -i WT_1_Aligned.sortedByCoord.out.bam -r /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/sacCer3.bed`

Output is an xls file and a summary txt file (mean & median values across all genes in sample).

Visualise TIN in boxplots in [Rstudio](https://github.com/friedue/course_RNA-seq2015/blob/master/03_mRIN.R) using ggplot

__Quality of RNA Seq Toolset (QoRTs)__
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

## Read Quantification
__Gene-based read counting__
- To compare the expression rates of individual genes between samples you need to **quantify the number of reads per gene.**
- Essentially you are **counting the number of overlapping reads**
- Need to clarify:
	- Overlap size (full read vs partial overlap)
	- Multi-mapping reads
	- Reads overlapping multiple genomic features of the same kind
	- Reads overlapping introns
![enter image description here](http://htseq.readthedocs.io/en/release_0.10.0/_images/count_modes.png)

**Tools to count reads:**

`htseq-count` has 3 modes union, intersection strict, and intersection nonempty (image above). Union mode is recommended which counts overlaps even if a read only shares part of its sequence with a gene but disregards  reads that overlap more than 1 gene. http://htseq.readthedocs.io/en/release_0.10.0/index.html

`featureCounts` counts reads if any overlap is found with a gene. Can exclude multi-overlap reads or include then for each gene that is overlapped. This is a package of Subread so need to `ml Subread`

`QoRTs` also does counting - Nobby uses this.

The underlying gene models supplied to the quantification program via GTF or BED files will affect the gene expression quantification.

Count reads per gene:
`featureCounts -a /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/sacCer3.gtf -o featureCounts_results.txt alignment/*bam`

2 output files:
featureCounts_results.txt has actual read counts per gene
featureCounts_results.txt.sumary gives quick overview of how many reads were assigned to genes. 

Count reads overlapping with individual exons:
`featureCounts -a /home/camp/ziffo/working/oliver/projects/rna_seq_worksheet/sacCer3.gtf -f -t exon -O -o featCounts_exons.txt alignment/*bam`

N.B. if the exon is part of multiple isoforms in the annotation file, featureCounts will return the read counts for the same exon multiple times. *n = number of transcripts with that exon*. Remove the multiple entries in the result file before going on to differential expression analysis.

**Preparing an annotation:**
To **assess differential expression of exons**, create an annotation file where overlapping exons of different isoforms are split before running featureCounts. Use `dexseq_prepare_annotation.py` script of DEXSeq package or `QoRTs`.

__Isoform counting__
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

# Normalising and Log Transforming Read Counts
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

**Tools for DGE:**
`edgeR` (better for false positives, less conservative, recommended if <12 replicates)
`DESeq`
`DESeq2` sample wise size factor
`limma-voom`
`cuffdiff`slow, cant support multifactored experiments, can detect differential isoforms, high false positives

![enter image description here](https://lh3.googleusercontent.com/LVvCl3GXhNzUx5lyTrHsr0z_ZmI0nb51TBiY1-53VifMuYW8HR9-X54sfLwoH5gFyqahHOm8_QaWhg "Comparison of DGE programs")

## Running DGE analysis tools in R

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

## Ontologies

**[Sequence Ontology](http://www.sequenceontology.org/browser/obob.cgi)**
There are >2,400 terms associated with sequences in the genome. Sequence Ontology defines sequence features used in biological annotations.
To search a defined sequence term use the [Sequence Ontology Browser](http://www.sequenceontology.org/browser/obob.cgi)
To quickly search Sequence Ontology, use grep on the raw data:
`URL=https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so-simple.obo`
`wget $URL`
`cat so-simple.obo | grep 'name: gene$' -B 1 _A 6` where $ represents the end of line; -B -A prints lines before & after a match.

**[Gene Ontology](http://geneontology.org/)**
Connects each gene to one or more functions.
3 sub-ontologies for each gene product:
- Cellular Component (CC): cellular location where product exhibits its effect
- Molecular function (MF): How does gene work?
- Biological Process (BP): What is the gene product purpose?
Searching GO: use http://geneontology.org/ or https://www.ebi.ac.uk/QuickGO/

GO Download page:
http://geneontology.org/page/download-annotations
- 2 files to download: 
1. definition (term) file
`wget http://purl.obolibrary.org/obo/go.obo`

2. association file
 `wget http://geneontology.org/gene-associations/goa_human.gaf.gz`. In GAF compressed format defined at http://geneontology.org/page/go-annotation-file-gaf-format-21
Contains both gene and protein IDs.

To search the function of a gene use the [GeneCards](http://www.genecards.org/) database to easily locate the gene by name.

### Gene Set Enrichment Analysis
- Identify common characteristics within a list of genes
- When using GO terms, this is called "functional enrichment"
- Most common variant is the ORA (over-representation analysis): 
	- examines genes in a list > 
	- summarises the GO annotations for each gene > 
	- determines if any annotations are statistically over-represented.

Multiple **GO enrichment tools** exist, many are difficult to use and outdated. The best are:
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
-   For doing this you can use the gene-level count table obtained from Kallisto.
