


> # Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/)   is a user-friendly algo which extract both gene and transcript level gene expression directly from fastq files using the raw fastq files. Then use  [Sleuth](https://pachterlab.github.io/sleuth/)  (also developed by Pachter lab) to perform differential gene and transcript expression analysis.

**SVD (singular value decomposition) analysis**

-   For doing this you can use the gene-level count table obtained from Kallisto. I wrote everything in R and I can send you some litterature which explains a bit the underlying math and idea. Also happy to speak about it over skype.

# Rapid Approach: Kallisto - Sleuth pipeline
Author = [Lior Patcher](https://en.wikipedia.org/wiki/Lior_Pachter)

[Kallisto](https://www.youtube.com/watch?v=94wphB3GKBM) quantifies transcript abundances. It pseudoaligns reads against a transcriptome (not genome). Simple count-based approaches underperform when determining transcript level counts as they disregard reads that overlap with more than one gene. If the genomic feature becomes a transcript rather than a gene it keeps many reads that would have been discarded.

Method of pseudoalignment: for each read the program aims to identify the target that it originates from using k-mers. By ignoring exactly where in the genome a read originates from it is much faster than normal alignment. This approach does not generate a BAM file (alignment file) but instead produce a measure of how many reads indicate the presence of each transcript. These use a **deBruikin graph** to assign reads to an isoform if they are compatible with that transcript structure.
![enter image description here](https://www.frontiersin.org/files/Articles/169488/fgene-06-00361-r2/image_m/fgene-06-00361-g002.jpg)

Schema of a simple deBruijn graph-based transcript assembly:
- Read sequences are split into all subsequence k-mers (here: of length 5) from the reads.
- A deBruijn graph is constructed using unique k-mers as the nodes and overlapping k-mers connected by edges (a k-mer shifted by one base overlaps another k-mer by k􀀀1 bases).
- The transcripts are assembled by traversing the two paths in the graph

Although much faster than alignment-counting routines it **cant detect novel isoforms**. However, instead of direct isoform quantification, you can glean more accurate answers from alternative approaches, e.g., quantification of exons (Anders et al., 2012) or estimates of alternative splicing events such as exon skipping, intron retention etc. (e.g., MISO [Katz et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037023/), rMATS [Shen et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280593/)).

The main limitations to assigning reads to transcripts are:
- annotation transcripts are inconsistent
- many isoforms with very different lengths
- anti-sense and overlapping transcripts of different genes

**Tools**:
Sailfish and more updated version Salmon
[Kallisto](https://www.nature.com/articles/nbt.3519)

## Kallisto Workflow
ml kallisto

Two steps:
1. Build Index (10mins to run)
2. Quantify Reads (10mins to run)

For an individual single or paired end fastq file you can run:
```bash
#set changable elements
##set reference transcriptome (cDNA)
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.fa
IDX=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.cdna.fa.idx
SAMPLE=CTRL_3
R1=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_depleted/${SAMPLE}.fq
##if you have paired-end data then set the 2nd read files for input
R2=PATH_TO_FASTQ_reverse_strand.fq

#build kallisto index
sbatch -N 1 -c 8 --mem 40 --wrap="kallisto index -i $IDX $REF"

#create directory for kallisto data 
mkdir -p kallisto
#set subdirectory for output called out
OUTDIR=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto
#run kallisto quantification with quant command. -o sets output directory. -b specifies the bootstap sample number.
##paired-end mode
kallisto quant -i $IDX -o $OUTDIR -b 100 $R1 $R2
##single-end mode (set fragment mean length & standard deviation - Illumina generates fragments 180-200bp - acurately determine this from a library quantification with an instrument such as an Agilent Bioanalyzer)
kallisto quant -i $IDX -o $SAMPLE -b 100 --single -l 187 -s 70 $R1
```
However for multiple fastq files use a For Loop:
```bash
# Create output folder
mkdir -p kallisto

# Exit this script on any error.
set -euo pipefail

# Set Reference transcriptome (not genome).
REF=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.fa
# Set index to build
IDX=/home/camp/ziffo/working/oliver/genomes/sequences/human/gencode.v29.transcripts.cdna.fa.idx

# Build kallisto index
sbatch -N 1 -c 8 --mem=40GB --wrap="kallisto index -i $IDX  $REF"

for SAMPLE in VCP CTRL;
do
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files.
        R1=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_depleted/${SAMPLE}_${REPLICATE}.fq
        # The kallisto output directory.
        OUTDIR=/home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/${SAMPLE}_${REPLICATE}
        # Run kallisto quantification in single-end mode.
        echo "*** Running kallisto on single end sample: $SAMPLE"
        kallisto quant --single -l 187 -s 70 -i $IDX -o $OUTDIR -b 100 $R1

        # Copy the abundance file to a proper name - i.e. remove the long path name so that it only contains information on the sample e.g. VCP_3.tsv
        cp $OUTDIR/abundance.tsv $OUTDIR.counts.tsv
    done
done

```

### Kallisto Output

Output files:
- abundance.txt
- abundance.h5 (large scale format form of txt file)
- run_info.json

The main output of Kallisto is the **abundance.tsv** file with columns:
```
target_id   length eff_length est_counts 	tpm
ERCC-00002  1061   891.059      18946      243099
ERCC-00003  1023   853.059      1452       19460.7
ERCC-00004  523    353.059      455        14734.5
ERCC-00009  984	   814.059      319        4480.29
```
Column 3 = eff_length = scales the transcript length by fragment length distribution . (transcript length - mean fragment length + 1)
Column 4 = est_counts = transcript abundance count
Column 5 = tpm = Transcripts Per Million

# Sleuth
https://github.com/pachterlab/sleuth
vignette('intro', package = 'sleuth')
help(package = 'sleuth')
https://pachterlab.github.io/sleuth/walkthroughs
https://pachterlab.github.io/sleuth_walkthroughs/pval_agg/analysis.html
https://groups.google.com/forum/#!forum/kallisto-sleuth-users

Sleuth runs differential analysis on kallisto counts (based on pseudo-alignment using alignment to the transcriptome index). Kallisto output includes transcript abundance estimates along with bootstraps.  Sleuth differential analysis model uses a general linear model, which accounts for biological & technical noise.

analyze the RNA-Seq dataset in order to obtain both gene-level and transcript-level differential expression results that are consistent with each other. 

Run Rstudio markdown pipeline:
`Sleuth.Rmd` located in `/Volumes/lab-luscomben/working/oliver/projects/vcp_fractionation/expression/sleuth`

The Sleuth object can be very large and needs to be run as an Rscript through CAMP cluster batch job:
```bash
ml Anaconda2
ml sleuth
source activate rtest
sbatch -N 2 -c 10 --mem 120G -t 12:00:00 --wrap="Rscript ~/working/oliver/projects/vcp_fractionation/expression/sleuth/sleuth_camp.R"
```
Can submit job as a text script:
```bash
#!/bin/bash -l
#PBS -N JOBNAME
#PBS -l walltime=5:00:00
#PBS -l nodes=2:ppn=10
#PBS -l pmem=6gb
#PBS -j oe

ml Anaconda2
ml sleuth
source activate rtest

cd ~/working/oliver/projects/vcp_fractionation/expression/sleuth

R --file=~/working/oliver/projects/vcp_fractionation/expression/sleuth/sleuth_camp.R
```

Debugging conda:
```bash
conda clean --all
conda update --all
conda update -n base -c defaults conda
conda install anaconda
conda install --channel https://conda.anaconda.org/bioconda r-sleuth
ml sleuth/0.28.0-foss-2016b-R-3.3.1
```

# Kallisto > DESeq2

Load Gene Annotaiton & Import sample details as for Sleuth

## Tximport 
[https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

Import transcript abundance estimates for samples using the tximport. rows = gene ID, columns = sample ID. 

kallisto abundance.h5 files can be imported by setting type to "kallisto". 


## Manipulate Count Files
Concatenate all counts into 1 file that can be used for DESeq DE analysis.
```bash
# In the kallisto output directory match sample names to conditions
CTRL=(CTRL_1 CTRL_2 CTRL_3)
VCP=(VCP_1 VCP_2 VCP_3)

#set up variables to help create count files
CTRL_FILES=''
for NAME in ${CTRL[@]}; do
    CTRL_FILES+="${NAME}/abundance.tsv "
done

VCP_FILES=''
for NAME in ${VCP[@]}; do
    VCP_FILES+="${NAME}/abundance.tsv "
done

#check names are correct:
ls $CTRL_FILES $VCP_FILES

#create the counts.txt file using all samples abundance.tsv files & remove unnessary columns
paste $CTRL_FILES $VCP_FILES | cut -f 1,4,9,14,19,24,29,34,39 > counts.txt

#Alternatively:
paste /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/VCP_*.tsv  /home/camp/ziffo/working/oliver/projects/airals/expression/D7_samples/kallisto/CTRL_*.tsv | cut -f 1,4,9,14,19,24,29 > counts.txt
```
This produces a Counts file that can be used as an input into DESeq. This file looks like:
```
target_id          est_counts  est_counts est_counts  est_counts  est_counts  est_counts
ENST00000472225.6    6.11331    1.28737     10.0567     2.23768     2.1543      0
ENST00000496578.3    3.12574    3.97029     8.18468     0           1.6547      1.49019
ENST00000624726.1    6.73189    3.0668      1.16443     2.39545     4.31444     1.26448
```
You can change the header to include the sample names.


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTI1OTU3MTk2MiwtNDQ2NTY2MTA2LDE3MD
cxMjgwMjMsNzU2ODE4ODY0LC0xNTcyOTc0OTA2LDE1MzM0MTA0
MTgsNjMxNjYyMl19
-->