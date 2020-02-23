
># Sequencing Data

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

### NCBI & Entrez
- ginormous data store of sequence data.
- Entrez is NCBIs primary text search integrating PubMed with 39 other databases
- Entrez enables you to access NCBI databases via the command line. 
- Accession number applies to the complete database record for that entity. Updates and revisions remain under the same accession number but with a different version number.
- Prior to 2015 NCBI versions began with GI.
- Entrez web API allows us to query NCBI data using the URL construct:  https://service.nih.gov/?param1=value1&param2=value2&param3=value3
- place `\` before each `&` to escape the command line `&` meaning or alternatively place ' ' around the whole URL.
- For example to access AF086833.2 in FASTA format via the command line: 

`curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=AF086833.2&db=nuccore&rettype=fasta'`

### Entrez Direct Tools
- Entrez direct tool simplifies the command line access e.g. `efetch`, `esearch`, `xtract`

```bash
efetch -db=nuccore -format=gb -id=AF086833 | head

#e.g. accession number AF086833 in Genbank format
efetch -db=nuccore -format=gb -id=AF086833 > AF086833.gb
#e.g. accession number AF086833 in Fasta format. 
efetch -db=nuccore -format=fasta -id=AF086833 > AF086833.fa
#efetch can take additional parameters and select a section of the sequence
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=3

#It can even produce the sequence from reverse strands
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=1
>gb|AF086833.2|:1-5 Ebola virus - Mayinga, Zaire, 1976, complete genome CGGAC

`efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=2`
>gb|AF086833.2|:c5-1 Ebola virus - Mayinga, Zaire, 1976, complete genome GTCCG
```

Note strand 2 is the reverse complement (not simple complement) represented by `c5-1` tag.

- searching Entrez Direct:
	- use the project accession number from the published paper to search the data that comes with it:
 `esearch -help`
 `esearch -db nucleotide -query accession_number`

- To fetch the data from a search pipe it to efetch: 
`esearch -db nucleotide -query accession_number | efectch -format=fasta > genome.fa` 
`esearch -db sra -query PRJNA*** | efetch -format runinfo > sequencing_data.csv`

- Manipulate the sequencing library if required:
isolate specific run IDs 
isolate single-end reads: `cat sequencing_data.csv | cut -f 1,16 -d , | grep SRR.*SINGLE | cut -f 1 -d "," > single_ids.tx`


## Downloading from Sequence Database Repositories
ml ncbi-vdb
ml fastq-tools
ml seqtk

### Accession numbers

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
5. On the [SRA page](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=384592) show 200 results & then click on Send Results to Run Selector. This creates the SRA Run Selector table which has useful information on each sequencing file.
	- [SRA Handbook](https://www.ncbi.nlm.nih.gov/books/NBK47528/)
	- SRA experiment IDs start with `SRS****`. Unique sequencing library for a sample.
	- SRA run IDs start with  `SRR****` and `ERR****`. Data file linked to the sequencing library.
6. On the [SRA Run Selector page](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_2827488_130.14.18.97_5555_1539957447_3051893689_0MetA0_S_HStore&query_key=3) click Accession List icon under Download. This downloads a text file of each Sequencing run ID.

![enter image description here](https://lh3.googleusercontent.com/MtirNCqJgzENksukOhrIZnqNwDpcXbbx_LTHMz17FZpe-p3Mn3iN0bsRAxsKgVbHJn3PjIOmyOOe5Q)

### fastq-dump

Now that we have the accessions, we can get the sequence files in fastq format using fastq-dump from the SRA toolkit

In terminal download each of the `SRR****` IDs in the txt file using the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) tool `fastq-dump` - this converts the SRA sequences to a FASTQ file.  In command line type `fastq-dump SRR1553607`. Generates `SRR1553607.fastq` file

The `while loop` allows you to do this for all SRR IDs in a single command and `sbatch` speeds up the process by parallelising the request: 
```bash
while read line; 
do 
	sbatch -N 1 -c 1 --mem 32 --wrap="fastq-dump --split-files --accession $line"; 
done < SRR_Acc_List.txt
```

Alternatively using the `ids.txt` file list of SRR run IDs & invoke `fastq-dump` on each ID:
```bash
fastq-fump -X 10000 --split-files SRR******` 

#This can be done in one go using
## make reads directory
mkdir -p reads
## run fastq-dump on each line of text file with xargs
cat ids.txt | xargs -n 1 fastq-dump -X 10000 --split-files reads
```

- **Paired-end reads** are concatenated by SRA (because that is how the instrument measures them). Therefore need to separate these into 2 different FASTQ files: 1 forward read; 1 backward read using the command: 
`fastq-dump --split-files SRR1553607`. Generates files `SRR1553607_1.fastq` & `SRR1553607_2.fastq`
- know the origin of each read (forward vs reverse) - encoded in read name - some analysis tools require combining the 2 files into 1. The **forward** read is **filename_1** and **backward** read is **filename_2**
- Need to process the read as Split Reads/files

**Using accession code to download fastq files directly**
- Copy the URL Address for the relevant FASTQ files from the column in the [SRA run selector table](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_110487706_130.14.18.97_5555_1538386558_3578612985_0MetA0_S_HStore&query_key=2). 
- use `wget` or `curl` tools to download e.g. `curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR197/SRR1972739/SRR
1972739.sra` or  	`wget link_copied_from the_website`
- These generate SRA files which need to be converted to FASTQ using: `fastq-dump -X 10000 --split-files SRR1972739.sra` we specified to download the first 10,000 reads using `-X 10000` and split the "spots" using `--split-files`. Spots are reads.
- If there are many samples then download summary (right click on TEXT) & copy link location: then in command line run: `wget -O samples_at_ENA .txt "<paste LINK >"` #the quotation marks are crucial 
- use 11th column of TEXT file (Fastq file top) to feed the URLs of different samples `cut -f11 samples_at_ENA . txt | xargs wget`

## Interrogate the downloaded data fastq files
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

## Edit fastq file names

It is very error prone to generate filenames from meaningless numbers e.g. SRR1553607. Better to rename into something meaningful e.g. VCP_D7_1. If there are many samples then use programming constructs to create filenames like so:

```
# Iterate over arrays and create file names from them.
CTRL=(SRR3191542 SRR3191543 SRR3194428)
for NAME in ${CTRL[@]}; do
    CTRL_FQ+="PATH/${NAME}.fq"
done

VCP=(SRR3191544 SRR3191545 SRR3194430)
for NAME in ${ZIKV[@]}; do
    VCP_FQ+="PATH/${NAME}.fq"
done
```

## Seqkit report on fastq
Installed [SeqKit](https://github.com/shenwei356/seqkit) & [csvtk](https://github.com/shenwei356/csvtk) and placed in PATH

Use seqkit to generate a report on fastq files:
For whole folder:`seqkit stat fastq_files/*` 
For individuals fastq files`seqkit stat *.fastq`

The report tells you:
- number of reads
- read lengths: min, avg, max
```
file       format  type    num_seqs      sum_len  min_len  avg_len  max_len
CTRL_1.fq  FASTQ   DNA    8,038,427  474,081,305       20       59       60
CTRL_2.fq  FASTQ   DNA    7,822,813  461,011,169       20     58.9       60
CTRL_3.fq  FASTQ   DNA    9,226,252  542,706,850       20     58.8       60
VCP_1.fq   FASTQ   DNA   12,475,620  736,312,401       20       59       60
VCP_2.fq   FASTQ   DNA   16,605,828  980,025,826       20       59       60
VCP_3.fq   FASTQ   DNA   13,406,845  789,756,869       20     58.9       60
```
 
 ## Optimise and manipulate fastq files (optional)
 
Convert FASTQ file into 3-column tabular format (1st col name/ID, 2nd col sequence, 3rd col quality) then add optional columns e.g. sequence length, GC content: `seqkit fx2tab --name --only-id --gc *.gz`

Extract a subset of sequences from a FASTQ file: ` seqkit sample --proportion 0.001 *.gz | seqkit seq --name --only-id > id.txt`
ID list file: `head id.txt` prints the relevant run IDs with specified proportion 0.1%
Search by ID list file: `seqkit grep --pattern-file id.txt *.gz > subset.fq.gz`
 
 Find degenerate bases in FASTQ file:  `seqkit fx2tab` converts FASTQ to tabular format & outputs the sequence in a new column.  `seqkit fx2tab --name --only-id --alphabet *.gz | csvtk --no-header-row --tabs grep --fields 4 --use-regexp --ignore-case --pattern "[^ACGT]" | c svtk -H -t cut -f 1 > id2.txt`
Use grep search to filter the table: `seqkit grep --pattern-file id2.txt --invert-match *.gz > clean.fa`
Or locate degenerate bases K & N: `seqkit grep --pattern-file id2.txt viral.1.1.genomic.fna.gz | seqkit locate --ignore-case --only-positive-strand --pattern K+ --pattern N+`

Remove duplicated sequences: `seqkit rmdup --by-seq --ignore-case *.fq.gz > uniq.fq.gz`

Locate a specific motif in FASTQ format: `seqkit locate --degenerate --ignore-case --pattern-file *.gz` use the degenerate flag to identify pattern of interest
 
## Fastq file names

### Create bash script to change file names:
create 2 text files: list of old_names & list of new_names
combine into 1 text file in 2 columns
```bash
paste old_names.txt new_names.txt > rename.txt
paste oldnames.txt newnames.txt | column -s $'\t' -t > rename.txt
sed 's/^/mv /' rename.txt | column -s $'\t' -t > rename.sh
```
start file with: `#!/bin/bash`
each line is a new `mv` command
save text file on CAMP but with ending as `.sh`
run script: `bash rename.sh`
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0NzgzOTUxNjAsLTIwMTI4NjgzNjEsLT
Y5OTIwOTQ1NywxNjY2NDEyMjMzLC0xMzM0OTA0NTQ3XX0=
-->