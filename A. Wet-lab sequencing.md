>  # Wet-lab RNA Sequencing

# Wet Lab phases
1. RNA isolation
2. Library Preparation
3. high-throughput sequencing

 # RNA extraction
 
* silica gel based membranes or liquid-liquid extractions with acidic phenol chloroform
* remove DNA and proteins. Improve with DNase.
* quality control: Aligent bioanalyser creates an **RNA integrity number (RIN)** is objective way of assessing RNA quality & degradation. 10 = intact; 1 = degraded. RIN of 8 is generally accepted threshold before proceeding to RNA seq. Uses elecetophoresis and looks for densitometry spike at 28S and 18S rRNA bands - ratio of 28S/18S = RIN.
 ![enter image description here](http://tlcr.amegroups.com/article/viewFile/286/596/2055)

# Library Preparation
* cDNA fragments 150-300bp —> hybridisation to flowcell (50-150 bp)
* small transcripts <150bp is lost in standard RNA-seq preparation
* mRNA enrichment: remove rRNA and tRNA by selecting polyA tails using oligodT beads OR removing rRNA with complementary sequences (ribo-minus approach —> retains unspliced RNAs).
 
# Stand sequencing
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

## Single read vs. paired end reads
* single read = determines the DNA sequence of just one end of each DNA fragment
* paired end = sequence both ends of each DNA fragment making pairing and directionality information available. 
	* Disadvantages: (i) 20% more expensive (ii) measures same fragment twice, thus at same genomic coverage will utilise only half as many unique fragments as single end sequencing. 
* for detecting de novo transcriptome assembly in humans need 100-200 x10^6 paired end reads.
![enter image description here](https://www.yourgenome.org/sites/default/files/images/illustrations/bioinformatics_single-end_pair-end_reads_yourgenome.png)

## Naming Sample Files

- make each attribute of the data be represented by a single, isolated region of the filename
- Utilise heirarchy of the samples. Start with the most generic information in the file name, then become more specific. E.g sample_replicate number_paired file (1 or 2)
- Dont mix sample information with replicate information
- keep as simple as possible


# Experimental Design
 
**Variability** in results:
* need **replicates** to capture breadth of isolate noise
* Technical replicates = repeat library preparations from the same RNA sample —> avoid batch effects & lane effects. Should multiplex same sample over different lanes of same flowcell.
* Biological replicates = parallel measurements on different samples i.e. RNA from independent cells/tissues. Most RNA-seq have 3 biological replicates but ideally need 6 per condition to improve statistical power.
 
**Artificial RNA spike-in control**
* used to accurately quantify absolute transcript concentration. 
* RNA of known quantities (92 transcripts) is used for **calibration** eg ERCC ExFold RNA Spike-In Control Mix. R package `erccdashboard`. Different spike-in controls are needed for each RNA type.
* dont use spike ins to normalise between different samples (they dont account for differences in amount of starting material).
 
**Blocking and randomise**
* Randomly choose which samples to treat and sample
* Block samples into groups based on known sources of variation (sex, weight, cell cycle status) - subexperiments in each block increases sensitivity.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE1NzU3OTc2OCwtMTY3NzQyMjU4NV19
-->