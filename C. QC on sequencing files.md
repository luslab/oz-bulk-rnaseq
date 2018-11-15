> # Quality Control (QC) on Raw Reads
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
`ml MultiQC`

make a folder to store the results: `mkdir fastqc_results`
run FastQC on each sequencing file: `fastqc SRR5* -o fastqc_results/`
if there are many fastq files needing FastQC then speed up by running as `sbatch`: `sbatch -N 1 -c 8 --mem 32 --wrap="fastqc SRR5* -o fastqc_results/"`

Look at the results: `ls fastqc_results/ERR458493_fastqc/`
Summarise multiple FastQC outputs using the [MultiQC tool](http://multiqc.info/):
[run `multiqc`](http://multiqc.info/docs/#running-multiqc) within the `fastqc_results` folder
Go to the folder with the fastqc files in and simply run: `multiqc .`
Open the MultiQC html report: `open multiqc_report.html` or [open in the browser](https://multiqc.info/docs/#using-multiqc-reports)
[Interpret](https://www.youtube.com/watch?v=qPbIlO_KWN0) the MultiQC report: 
general stats
quality scores (Phred score)
Over-represented sequences: these indicate how much % the adapter sequences compose

### Filter low quality bases & Trim adapters (Essential Step)
There are many QC tools available (most in bash but some in R - bioconductor) each with basic QC methods plus unique functionality. Best ones include:
`Trimmomatic`  application note in Nucleic Acid Research, 2012, web server issue
`BBDuk` part of the BBMap package
`FlexBar` Flexible barcode and adapter removal  published in Biology, 2012
`CutAdapt`  application note in Embnet Journal, 2011 - advised by Nobby
`Trim Galore` is a wrapper around `Cutadapt` - this performs quality trimming then adaptor trimming all in 1. It can also run FastQC.

`ml Trim_Galore`
`ml cutadapt`
`ml trimomatic`
`ml FASTX-Toolkit`
`ml Bowtie2`

1. `mkdir trimmed_results`

2. Run Trim Galore 
`trim_galore -q 20 --length 20 --gzip -o /home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_results /home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/SRR5483788_1.fastq`

for multiple sequences can parallelise by using for loop & sbatch:
```bash
for file in ~/working/oliver/projects/airals/fastq_files/D7_samples/SRR5*_1.fastq
do
	sbatch -N 1 -c 1 --mem 32 --wrap="trim_galore -q 20 --length 20 --gzip -fastqc -o /home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_results $file";
done
```

`-q 20` = trim low quality ends Phred <20 & remove adapters
`--length 20` = discard reads shorter than 20bp
`-o` = specify output directory
`-fastqc` = will perform FastQC after low quality & adapter removal
You can specify the precise adapter sequence you want removed e.g. `-a AGCGCTAG` but if this isnt specified explicitly, Trim Galore auto-detects the Illumina adapter sequence that was used. It specifies what it removed in the output txt file.

@Raphaelle used:
adapter removal: `fastx_clipper -Q 33 -l 24 -a $adapt -i ${paths}${file} > ${paths}${data}-clipped.fastq`

3. Re-run MultiQC
run `multiqc` within the `trimmed_fastqc_results` folder
Go to the folder with the trimmed fastqc files in and simply run: `multiqc .`
Compare this new Trimmed MultiQC HTML report with the report on the Raw FastQC prior to trimming.

### Remove rRNA & tRNA (Optional step)
Can assess after mapping to see how much RNA has mapped to ribosomal RNA reference genome. If >5% then consider depleting these reads. 
There are disadvantages to removing these sequences:
Generally for most RNA-seq remove rRNA at this point pre-processing.

1. Create rRNA&tRNA reference genome
`ml BEDOPS`
```bash
awk '/rRNA|tRNA|Mt_rRNA|Mt_tRNA|miRNA|misc_RNA|snRNA|snoRNA/{print $0}' gencode.v28.primary_assembly.annotation.gtf > temp.gtf
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' temp.gtf > gencode.v28_ribosomal.gtf
gtf2bed < gencode.v28_ribosomal.gtf > gencode.v28_ribosomal.bed
sed -i -r 's/^([0-9]+|[XY]|MT)/chr\1/' gencode.v28_ribosomal.bed
```

2. Index the rRNA reference genome

ml Bowtie2
ml BEDTools
ml SAMtools
ml BWA

```bash
mkdir rRNA_depleted

#BEDTools (any many other programs) need TABS, not spaces.
#sed 's/  */\t/g' Homo_sapiens.GRCh38.77_ribosomal.bed > temp.bed
#mv temp.bed Homo_sapiens.GRCh38.77_ribosomal.bed
bedtools getfasta -fi /home/camp/ziffo/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa -bed gencode.v28_ribosomal.bed -fo gencode.v28_ribosomal.fa
bwa index gencode.v28_ribosomal.bed
bowtie2-build gencode.v28_ribosomal.fa gencode.v28_ribosomal
samtools faidx gencode.v28_ribosomal.fa
```

3. Map sequences to the rRNA & tRNA genome

```bash
#set shortcuts
##set the desired name of the output fastq file with the RNA aligned reads depleted
OUTPUT=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/rRNA_depleted/SRR5483788_1_no_rRNA.fq
##set reference genome as the ribosomal genome
REF=/home/camp/ziffo/working/oliver/genomes/annotation/GRCh38.p12/gencode.v28_ribosomal
##set the input fastq file
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/trimmed_results/SRR5483788_1_trimmed.fq.gz
##set name of output SAM file with alignment info to the RNA genome
RNA_SAM=/home/camp/ziffo/working/oliver/projects/airals/fastq_files/D7_samples/rRNA_depleted/SRR5483788_1_with_rRNA.sam

#map to ribosomal genome
sbatch -N 1 -c 8 --mem 40 --wrap="bowtie2 -q -p 8 --un $OUTPUT -x $REF -U $FASTQ -S $RNA_SAM"
```
-q = input is fastq; -p 8 = launch 8 alignment threads; --un (path) = write unpaired reads that didnt align to this path; -x bt2 = index filename prefix; -U file.fq = files with unpaired reads (can be .gz); -S sam = sam output file

4. Re-run FastQC & MultiQC step
`ml pandoc`
```bash
for file in ~/working/oliver/projects/airals/fastq_files/D7_samples/rRNA_depleted/SRR5*_1.fastq
do
	sbatch -N 1 -c 8 --mem 40 --wrap="fastqc $file -o rRNA_depleted_fastqc_results/";
done
run multiqc within the rRNA_depleted_fastqc_results/ folder: multiqc -f .
```
 
### Error Correction (Optional step)
Nobby doesn't routinely do this. 
Sequencing instruments occasionally make random errors. Can recognise errors by computing k-mer density. K-mers are all possible short subsequences of length `k` contained in a string (sequence) that occur many times. They can be 2-5 bases long. K-mers that contain errors will be rare so can be used in error correction but also genome identification/classification & alignement.

- The 2 base long k-mers ( 2-mers ) are AT , TG , GC and CA
- The 3 base long k-mers ( 3-mers ) are ATG , TGC and GCA
- The 4 base long k-mers ( 4-mers ) are ATGC , TGCA
- The 5 base long k-mer ( 5-mer ) is ATGCA

FastQC error correction programs correct or remove reads that appear to have errors in. Very useful when data is of high coverage.

`ml BBMap`
`bbmap` package using `tadpole.sh` error corrector
`tadpole.sh in=SRR5*_1.fastq out=tadpole.fq mode=correct`
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTU5OTMyMjkwNCwtMTg5OTgyMDIyXX0=
-->