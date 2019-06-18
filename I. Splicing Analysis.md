


> # Differential Splicing Analysis
[Video RNA seqblog](https://www.rna-seqblog.com/differential-splicing-analysis-with-rna-seq-current-applications-approaches-limitations/)
[Coursera tools video](https://www.coursera.org/lecture/genomic-tools/tools-for-transcriptomics-2-rna-seq-ESKKW)
![enter image description here](https://lh3.googleusercontent.com/mAVrLXkGSEAroRb6hIXjSFkHJhrHrXo_T0wKWqVrTXofmVU1I_NTXU70mhs3LG6OuvMPcdzSLz9JyQ)


![enter image description here](https://lh3.googleusercontent.com/BGrCH8FUIUWYDSY2U1b_kBzc0AtCXguL7EQdOR3ozQ1wCLdZPvjk58ImA96UpRVRJWUX4-b5IjuyhA)

Cassette exon = Exon skipping
Alternative splice sites incorporate different regions of an exon in the mRNA isoform
Mutally exclusive exons = include one or the other exon
Concatenate = link exons into a chain to form mRNA

Alternative splicing is used normally to enable a gene to produce different proteins which maybe tissue specific. However aberrant splicing can cuase disease by inappropriate exon skipping (e.g. BRCA1 exon 18 in cancer) and intron retention.
![enter image description here](https://lh3.googleusercontent.com/jeQ92hUjro9rmqVMyY53dnH1cGWIRAN466rpSTAceqTWdWFm8xN0aO3DgUfebTQZsvWNWg8UhYEgXA)

In Differential Gene Expression you compare mutant vs control read counts for each gene. 
In **Differential Splicing Expression** you compare all the different mRNA isoforms between mutant and control - there are many isoforms for each gene so there are many comparisons.
![enter image description here](https://lh3.googleusercontent.com/VyUh9L2RvsCzGjGbh1fV5pPq2gGCU4gadhSaCk7BvpvceELVaUiMngWYwXIwFhbvU1xTQmAR0f4DFg)

In **Differential Isoform Usage** you compare the usage of all the alternatively spliced transcripts **RELATIVE** to the total gene expression between mutant & control. I.e. take the Isoform expression counts divide by that gene expression counts. This is the important calculation with Diffential Splicing Analysis.
![enter image description here](https://lh3.googleusercontent.com/yc1EU-oIRVA1EuDyUkmW0rQewMRlBPWaOtgLIkpdAxMqaMr74KSbzfc4kekQKSEpV62ZaWv0WX1ElA)
![enter image description here](https://lh3.googleusercontent.com/vCWCaScTZWxndCmsKzKNrrQULS51WUsJNo-e64d5YP-s1m6kfWd_0pyNPBsItwIlo-Rhaz1O7N5vpw)



# Approaches to Differential Splicing Analysis
**1. Isoform Resolution Approach**
Differences in complete isoform proportions (expression) between samples. Does NOT rely on transcriptome annotation as with the other approaches.
Method: Assemble all full length isoforms > Quantify expression of each isoform > Test differeneces in relative abundance
Advantages: investigate full length isoforms. Doesnt rely on transcriptome annotation > better for identifying novel events.
Limitations: complex, ambiguity with different reads aligning to same position.
Tools: **Cuffdiff**, DiffSplice, SplicingCompass

**2. Exon Usage Approach**
Differential exon-level expression between conditions. Compares individual exon expression count.
Method: Avoids complexity of transcript assembly & expression estimation. Assumes differential usage of non-terminal exons equates to differential splicing. **Requires transcriptome annotation** to identify all possible exon units for each gene. Genes models are flattened into exon counting bins (tourquoise in image below). Overalpping exons and exons with alternative boundaries are split into separate bins (dotted lines in image). Then quantifies individual exon bins. Then tests differential usage of each exon bin using generalised linear model (same as with differential gene expression)
Tools: **DEXSeq**, DSGSeq, GPSeq, SOLAS
Advantages: easier as doesnt resolve full length isoforms, doesnt make an abundance estimation at the transcript level
![enter image description here](https://lh3.googleusercontent.com/00opX631NuAn6nJrNQatEU2G9n6Hk-e0UxaMGqVwGV6vJUI3VHrUEaQ3CPcnd1DpIqpqEFpoVV8uaA)
![enter image description here](https://lh3.googleusercontent.com/yLcLWxmo7DlPLJyzSeLxqlae97F9a69sXdXJeDOxf2ct-_e7wj9iiNcWAxF6hMC4UQccwrSbR-r7gQ)

**3. Splicing Event Approach**
Differential inclusion of alternative splicing events between samples. Focus on localised alternative splicing patterns. Builds on exon-based approach by incorporating splice site information. Compresses  alternative transcripts into a single unit:
![enter image description here](https://lh3.googleusercontent.com/mxwnemNNIlCnwrJdpNKroGcigW-yZHFv0C8jIhXaSdxHhoX7eCElMdybH2mi7DDZySNQzmIx_fir0Q)

Method: Differential splicing is assumed from **differential inclusion** of a particular splicing event (5 simple AS events):
![enter image description here](https://lh3.googleusercontent.com/Yp1gfeFb5Vv45TQHCfigT4J29TwTpT4NqcS9QUNAPB2Zp9XvK-uE1ZLqpSv-m8znWoYIilXvpKu6RA)
Identifies the different alternative splicing patterns from transcriptome annotation. Exons from all transcripts with the same flanking exons and binary event exon are collapsed into single units (torquoise below).

Splicing event is quantified as **percent spliced in (PSI)**. 
PSI = Alternative Splicing Event / Normal Splicing Event. 0 < PSI < 1. Image below PSI 0.76 = 76% of the genes expression is coming from the exon inclusion event.
**PSI index = (a + b)/(a + b + 2c)**, where a and b = the number of splice-junction reads connecting the alternative exon to the upstream and downstream constitutive exons, respectively, and c = the number of junction reads connecting the two constitutive exons.
Calculates change in inclusion (delta PSI) for all potential splicing events between conditions. Tests if delta PSI exceeds a set threshold e.g. 5% (likelihood-ratio test). Visualise results in a sashimi plot.
**Advantages**: doesnt make assumptions that is required with full isoform resolution, 
**Limitations**: only tests the 5 simple classical splicing patterns (skipped exons, retained introns, alt 5' splice site, alt 3' splice site, mut excl splice sites )& **misses complex patterns not classified** (compound events eg Alt 5'ss + ES + Alt 3'ss, multiple skipped exons, non binary splicing - account for ~20% of AS events); splicing events are **strictly binary** (inclusion/exclusion of event):
![enter image description here](https://lh3.googleusercontent.com/IdaivGvJ6g2Y4oNTLGwnrY4pqhbv_wemYzT_sYYX71_fWnTn9jNp1bLZcLKM6WSuK0RVzaGiYDhR4Q)

![enter image description here](https://lh3.googleusercontent.com/D52SoO2TmOgf9tzlCamVAVRvDkZGGqoBAO1zJbcz-bNg6zntuJ2T8VUrUs7t6NtPDc4grPlkyKihZw)
Tools: **VAST-TOOLS**, rMATS, JuncBASE, JETTA, SpliceSeq

![enter image description here](https://lh3.googleusercontent.com/k05ELwI_1WON11e6TULfOIQ8Mc9w-248qVu2LkCJ148MrOcK3wA9EsVOqJcTWus24yGIVNV3zRugVg)

![enter image description here](https://lh3.googleusercontent.com/V1TmFNGiyJK7TehC_I6G0n_G9l71OA282Q1RXPD1qEJCLUC1xB8cFQWB1jIEJEpI_1sPk5uolLAGmA)

![enter image description here](https://lh3.googleusercontent.com/Y8GI4VPHpZ6V6Rp99jtcERChY83IWskw6PC-4Fc4NuIUOUlKbJJ-VLrky9bhcuEWWn7VNq31jloC7A)

**4. Alternative Splicing Graphs Approach**

Method: build splicing graph from transcriptome annotations, representing all possible AS variants. Then identify & quantify AS events.
Advantages: utilising read alignment information, can identify novel splicing events (not possible with DEU & AS event approaches). Can annotate & quantify complex splicing patterns.
Tools: ASGAL, Whippet, SGSeq
![enter image description here](https://lh3.googleusercontent.com/vLTaa0NevRDtwSSIIHXQusqYmSolMCtMYJPLFzc53zC3SfoKh_mnNwSs61W0iAla6trRXRFo5ifH2A)

# Limitations of Differential Splicing Analysis

## 1. Overlapping Transcriptome Features
Transcriptome features found at same gene location: completely overlap, partially overlap, shared boundaries. Introduces complexity & ambiguity for splicing analysis. 
DNA is double stranded: the 2 strands can produce different isoforms. But even on same strand you can produce different isoforms. Can work out which reads come from which strand:
![enter image description here](https://lh3.googleusercontent.com/sO5sJcsyqVa7u2FfuTooa5YOUbIM7X2Vali9zWR45BBUj0m3CIh4eOFCPLvL3dMSs0eFeaj7vWN5QA)

Deal with opposite strand issue with strand-specific library prep:

Isoform resolution approach deals with same strand isoforms calculate the **maximum likelihood estimation** to determine probability of a read belonging to a particular isoform. Overlapping genes are a major problem: they increase the number of transcripts & potential for overlapping exons.

Exon based approach deals with this by aggregating overlapping genes into a single group and testing differential exon usage within the group. Then removals all overlapping exons and tests remaining structures. This looses exons. For overlapping genes it has ambuguity at gene level for differential splicing - mistake differential gene expression for differential exon usage.

Splicing event approach deals with this by identifying splicing events based on transcriptome annotation. Then annotates events separately for each gene. Then tests each event of each gene independently (doesnt group genes together). The disadvantage of this is that ambiguious reads aligning to more than 1 gene will be counted towards the expression of both isoforms.

![enter image description here](https://lh3.googleusercontent.com/iHZ1u2hZF4WkZZW_D4AOviQO9rtSTJ_uemJ9fYQ_W3ibJ8xpdYOv0np1GDvcY7XNaONrc1ELabmhCA)

## 2. Complex Splicing Events

Neglected by the differential exon usage & alternative splicing event approaches.

![enter image description here](https://lh3.googleusercontent.com/IdaivGvJ6g2Y4oNTLGwnrY4pqhbv_wemYzT_sYYX71_fWnTn9jNp1bLZcLKM6WSuK0RVzaGiYDhR4Q)


# Tools
[http://www.rna-seqblog.com/tag/alternative-splicing/](http://www.rna-seqblog.com/tag/alternative-splicing/)

- [VAST-TOOLS](https://github.com/vastgroup/vast-tools): Ben Blancoe's lab. Used by Raphaelle. Looks at intron retention & 2 junction reads within each exon to look at Microexons.
- [Matt](https://academic.oup.com/bioinformatics/article/35/1/130/5053311): UNIX command line tool. Downstream analysis of VAST-Tools PSI output table to provide exon comparisons; motif RNA maps; [http://matt.crg.eu/](http://matt.crg.eu/)
- [PSI Sigma](https://github.com/wososa/PSI-Sigma): A new [PSI index](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848607/). Traditionally, the PSI index is denoted as (a + b)/(a + b + 2c), where a and b = the number of splice-junction reads connecting the alternative exon to the upstream and downstream constitutive exons, respectively, and c = the number of junction reads connecting the two constitutive exons. We modified the PSI index as follows: 
![enter image description here](https://lh3.googleusercontent.com/7AkyPb4-l3ftfBmhfl0gqHpvCB85STqY3iAAiFf7PGtVUD1x4KrcIYKC_NJNGQEVAoIka-aLPKVq6A)
C1 and C2 = the upstream and downstream constitutive exons, respectively. 
C1Si =the total number of junction reads whose 5′ splice site is connected to the upstream constitutive exon in a given splicing event.
C2Sj = the junction reads whose 3′ splice site is connected to the downstream constitutive exon.

- [SpliceDetector](https://www.nature.com/articles/s41598-018-23245-1): SpliceGraph forms based on freq. of active splice sites in pre-mRNA. Then, compares transcript exons to SpliceGraph exons. Discovers AS events from known transcripts. Simple & Fast. Transcript ID > build SpliceGraph using Exon coordinates > identify AS events
- [ASGAL](https://asgal.algolab.eu/): predicts events that use splice sites which are novel with respect to a splicing graph.  Directly align reads to a splicing graph.
- [Portcullis](https://github.com/TGAC/portcullis): removes invalid splice junctions from pre-aligned RNA seq data. Splice aware aligners often produce many false positive splice junctions. Filters culls splice sites which are unlikely to be genuine. 

- [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html): focused on differential exon usage. [Vignette](http://127.0.0.1:12657/library/DEXSeq/doc/DEXSeq.pdf).
- JunctionSeq is like DEXSeq with junction reads included (and is written by the QoRTs team). JunctionSeq vignette - they have a great walkthrough that ... walks you through the whole process from beginning to end inc. QoRTs
- [rMATS](http://rnaseq-mats.sourceforge.net/): useful for comparing with other ENCODE datasets. only lists novel events that use annotated splice sites
- [MISO](http://genes.mit.edu/burgelab/miso/):Bayesian inference to estimate the probability for a read to be issued from a particular isoform. supplies confidence intervals (CIs) for: (i) estimating of exon and isoform abundance, (ii) identifying differential expression. It can be applied for analyzing isoform regulation.
- [MAJIQ] is also good but parsing the output is a bit annoying (but the default was the best looking one of the lot!). quantifies the relative abundances of a set of Local Splicing Variations which implicitly represent combinations of AS events involving both annotated and novel splice sites
- Whippet is new and lightweight, but you can't really see what it is up to or the reads it has aligned (edited)
- [LeafCutter](https://www.nature.com/articles/s41588-017-0004-9) https://github.com/davidaknowles/leafcutter quantifies differential intron usage across samples, allowing the detection of novel introns which model complex splicing events. requires as input the spliced alignments. 
- [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html#overview-of-alternative-splicing-workflow)
- [IR Finder](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4)
- SUPPA2: only able to detect AS events that are in the annotation. requires the quantification of the input transcripts, which can be obtained by using Salmon
- [DARTS]([https://github.com/Xinglab/DARTS](https://github.com/Xinglab/DARTS) [https://www.nature.com/articles/s41592-019-0351-9](https://www.nature.com/articles/s41592-019-0351-9) uses deep learning to analyse alternative splicing
- [Cuffdiff](http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff), [ALEXA-seq](http://www.alexaplatform.org/alexa_seq/), [MISO](http://genes.mit.edu/burgelab/miso/), [SplicingCompass](http://www.ichip.de/software/SplicingCompass.html), [Flux Capacitor](http://flux.sammeth.net/capacitor.html), [JuncBASE](http://compbio.berkeley.edu/proj/juncbase/Home.html), [DEXSeq](http://bioconductor.org/packages/2.11/bioc/html/DEXSeq.html), [MATS](http://rnaseq-mats.sourceforge.net/), [SpliceR](http://www.bioconductor.org/packages/2.13/bioc/html/spliceR.html), [FineSplice](http://nar.oxfordjournals.org/content/early/2014/02/25/nar.gku166.full), [ARH-seq](http://nar.oxfordjournals.org/content/early/2014/06/11/nar.gku495.full),

**Alternative splicing, alternative expression**[[24447644](http://www.ncbi.nlm.nih.gov/pubmed/24447644),  [24885830](http://www.ncbi.nlm.nih.gov/pubmed/24885830),  [24058384](http://www.ncbi.nlm.nih.gov/pubmed/24058384),  [24549677](http://www.ncbi.nlm.nih.gov/pubmed/24549677),  [24951248](http://www.ncbi.nlm.nih.gov/pubmed/24951248),  [25511303](http://www.ncbi.nlm.nih.gov/pubmed/25511303)]

Cuffdiff [[23222703](http://www.ncbi.nlm.nih.gov/pubmed/23222703)], DEXSeq [[22722343](http://www.ncbi.nlm.nih.gov/pubmed/22722343)], ALEXA-seq [[20835245](http://www.ncbi.nlm.nih.gov/pubmed/20835245)], IUTA [[25283306](http://www.ncbi.nlm.nih.gov/pubmed/25283306)], FineSplice [[24574529](http://www.ncbi.nlm.nih.gov/pubmed/24574529)], PennSeq [[24362841](http://www.ncbi.nlm.nih.gov/pubmed/24362841)], FlipFlop [[24813214](http://www.ncbi.nlm.nih.gov/pubmed/24813214)], SNPlice [[25481010](http://www.ncbi.nlm.nih.gov/pubmed/25481010)], spliceR [[24655717](http://www.ncbi.nlm.nih.gov/pubmed/24655717)], GESS [[24447644](http://www.ncbi.nlm.nih.gov/pubmed/24447644)], RNASeq-MATS [[23872975](http://www.ncbi.nlm.nih.gov/pubmed/23872975)], SplicingCompass [[23449093](http://www.ncbi.nlm.nih.gov/pubmed/23449093)], DiffSplice [[23155066](http://www.ncbi.nlm.nih.gov/pubmed/23155066)], SigFuge [[25030904](http://www.ncbi.nlm.nih.gov/pubmed/25030904)], SUPPA [bioRXiv], CLASS [bioRXiv], SplAdder [bioRXiv], SplicePie [[25800735](http://www.ncbi.nlm.nih.gov/pubmed/25800735)].

Annotated versus novel exploratory events. 

## Analysis approach
https://www.nature.com/articles/nmeth.1503.pdf

Infer structural infromation about the transcript 
Infer the strand by examining splice site spanning reads
Each transcript isoform has very few exons & exon-exon junctions that are unique to that isoform.

Can ignore structure of full length transcript & focus on individual sequence features.

![enter image description here](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=info:doi/10.1371/journal.pcbi.1004393.g006)
yellow gene = noncoding RNA gene.
brown & green genes = coding genes
Few exon-exon spanning genes.

# Gene Isoform counting

Gene isoforms are mRNA produced from the same locus but with different protein coding sequences.  5 modes of alternative splicing are recognised:
1.  Exon skipping or cassette exon.
2.  Mutually exclusive exons.
3.  Alternative donor (5') site.
4.  Alternative acceptor (3') site.
5.  Intron retention.

# VAST-TOOLS

## Overview

VAST-tools output = any splicing changes OVER TIME. i.e not just retained introns
Focus is how splicing patterns changed over time (Day 0; Day 7; Day 14 & Day 21) in the CTRL & VCP groups. Rather than directly compare splicing differences at each stage between CTRL & VCP, I compared the groups of genes at each stage which were exhibiting changes in splicing over time in control group versus VCP group. 
Most changes in IR in CTRL occurred at day 14 whilst in VCP the same events were occurring at day 7, premature IR/splicing in VCP mutants.

The `Splicing_VASTOOLS.sh` script is located in `/home/camp/ziffo/working/oliver/scripts/intron_retention`

There are 6 steps:
1. Align (this is redone in VAST-TOOLS using bowtie - needs loading)
2. Combine outputs into 1 summary table
3. Merge technical repeats (optional)
4. Tidy (optional)
5. Differential Splicing analysis
6. Plot the output

## Alignment
ml R
ml Bowtie

https://github.com/vastgroup/vast-tools#alignment

Use untrimmed fastq files - use raw reads. Define reference genome species (Hsa = human). 
```bash
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools
# Not necessary to create output folder as VAST-tools auto creates an output folder within your current working directory named "vast_out" and within that is "tmp" and "to_combine" folders

#set shortcuts
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D0_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D7_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D14_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D21_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D35_samples/SRR*_1.fastq
FASTQ=/home/camp/ziffo/working/oliver/projects/airals/reads/D112_samples/SRR*_1.fastq

#run vast-tools on each FASTQ file separately at each time point. Dont specify output as all files need to be in same subfolder > output auto goes into a folder called vast_out. Run from the vast-tools directory
for SAMPLE in $FASTQ
do
	sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools align $SAMPLE -sp Hsa"
done
```

## Merging Outputs
https://github.com/vastgroup/vast-tools#merging-outputs

Merge the aligned output files for technical replicates when read coverage for independent replicates is not deep enough for a complete AS analysis.  Ideally have >150 million reads per sample for VAST-TOOLS AS analysis.  Raphaelle merged the 3 samples for each of VCP & CTRL at each time point.
If no technical replicates then skip this.

```bash
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools
#Prepare config_file from sample group txt file (needs 2 columns: fastq file name & group separated by a tab)
awk '{print $3"\t"$2}' /home/camp/ziffo/working/oliver/scripts/intron_retention/VASTOOLS_merge_groups.txt | tail -31 > /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/config_file

CONFILE=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/config_file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out

#Run it -- takes about 10 minutes on the 31 samples
vast-tools merge --groups ${CONFILE} --o $OUT --sp Hsa --move_to_PARTS
```

## Combining results
ml R
https://github.com/vastgroup/vast-tools#combining-results

About 2G of memory required; completed in about 10 min
Combine aligned files that are stored in the folders `to_combine`  to form one final table called `INCLUSION_LEVELS_FULL-Hsa6-hg19.tabz` (if using old legacy version) or 5 tables in v.2.0.0. This is the table that you send to differential splicing command. Can specify hg38. The output directory contains the sub-folders to combine. 
```bash
#move to to_combine directory with merged aligned files
cd /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out/to_combine
#set aligned output file
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out

#  create the old legacy version INCLUSION_TABLE.tab single output then specify `--noANNOT`
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa -v --noANNOT"

# run vast-tools combine using new v2.0.0 ANNOT tool - identifies annotated exon-exon reads
sbatch -N 1 -c 8 --mem=40GB --wrap="vast-tools combine -o $OUT -sp Hsa --IR_version 2"
# This produces 5 INCLUSION_TABLE files in the raw_incl folder. 
##ANNOT = identifies & profiles annotated exons. This is the novel feature of v.2.0.0
##COMBI = splice site based pipeline
##EXSK = only alternative Exons are considered. Transcript based pipeline, single.
##MIC = microexon pipeline
##MULTI = transcript based pipeline, multiexon

#check output
## old legacy table
cd vast_out
head INCLUSION_LEVELS_FULL-Hsa6-hg19.tab
## new v2.0.0
cd vast_out/raw_incl
head INCLUSION_LEVELS_ANNOT-Hsa14-n.tab
head INCLUSION_LEVELS_COMBI-Hsa14-n.tab
head INCLUSION_LEVELS_EXSK-Hsa14-n.tab
head INCLUSION_LEVELS_MIC-Hsa14-n.tab
head INCLUSION_LEVELS_MULTI-Hsa14-n.tab
```
Format of the combine output is as follows: 
https://github.com/vastgroup/vast-tools/blob/master/README.md#combine-output-format

## Compare Groups & Differential Splicing Analysis
https://github.com/vastgroup/vast-tools#comparing-psis-between-samples
https://github.com/vastgroup/vast-tools#differential-splicing-analysis

-   `compare`: pre-filters the events based on read coverage, imbalance and other features, and simply compares average and individual dPSIs. It looks for non-overlapping PSI distributions based on fixed dPSI cut-offs. For more than 3 replicates, it is likely to be too stringent.
-   `diff`: performs a statistical test to assess whether the PSI distributions of the two compared groups are significantly different. It is possible to pre-filter the events based on the minimum number of reads per sample, but subsequent filtering is highly recommended (e.g. overlapping the results with the output of  `tidy`). For more than 5 samples per group it may also be over stringent.

As CAMP R module doesn't have psiplots R package need to create a conda environment to install packages. Once this environment has been made it is saved for future and can be activated using `source activate`. Then load relevant R module tools within the conda environment.
```bash
ml Anaconda2
ml R/3.5.1-foss-2016b-bare

source activate rtest
R
library("ggplot2")
library("optparse")
library("MASS")
library("RColorBrewer")
library("reshape2")
library("grid")
library("parallel")
library("devtools")
library("psiplot")
#quit R in cluster
q()
#save workspace image
```
To create a conda environment from scratch:
```bash
conda create -n rtest r-essentials r-devtools
source activate rtest
# install the package normally by calling R
R
install.packages("optparse")
install.packages("ggplot2")
install.packages("MASS")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("devtools")
install.packages("psiplot")
devtools::install_github("kcha/psiplot")

#to deactivate environment
conda deactivate
```
run `vast-tools compare` and `vast-tools diff` in this rtest conda environment outside of R

There are 2 approaches to comparing samples: 
1. Time effect (delta change in differential splicing between different time points of iPSC differentiation in CTRLs and VCPs independently). Effectively is a pair-wise analysis like repeated measures but between 2 time points only.
2.  Mutant effect (VCP vs CTRL at specified time points)

### Time Effect

```bash
mkdir /home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/time_effect
INFILE=/home/camp/ziffo/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out/INCLUSION_LEVELS_FULL-Hsa14-hg19.tab
```
Run the Compare TIME-EFFECT section of the script in `/home/camp/ziffo/working/oliver/scripts/intron_retention/Splicing_VASTOOLS.sh`
This performs the VAST-tools Compare & Diff commands on each of the WT & VCP time-point comparisons

### Mutant Effect

Run the Compare MUTATION-EFFECT section of the script in `/home/camp/ziffo/working/oliver/scripts/intron_retention/Splicing_VASTOOLS.sh`

Can use VAST-TOOLS here to calculate differentially expressed genes: `compare_expr`
Output file is created in directory of input file. This reports the differentially spliced AS events between the 2 groups (based on difference in average inclusion levels - delta PSI)

Output file of diff command = INCLUSION-FILTERED.tab

### View diff output files
```bash
order tab file by MV value:
more INCLUSION-FILTERED.tab | sort -k6 -r | awk '{ if ($6 >= 0.2) { print } }' | awk '{ if ($5 >= 0) { print } }'

# count number of genes with MV > X:
awk '{ if ($6 >= 0.2) { print } }' INCLUSION-FILTERED.tab | wc -l
# count number of genes with +ve intron retention (VCP vs CTRL):
awk '{ if ($6 >= 0.2) { print } }' INCLUSION-FILTERED.tab | awk '{ if ($5 >= 0) { print } }' | wc -l
```


# Matt
[http://matt.crg.eu/](http://matt.crg.eu/)

Two phases: Data preparation > Data analysis.
Input files are tab separated plain text tables

The workflow is as follows: pre-process the output table of VAST-TOOLS for extracting PSI values with [get_vast](http://matt.crg.eu/#get_vast), define groups of events to be analyzed with [def_cats](http://matt.crg.eu/#def_cats), and eventually run a high-level analysis (50 intron features).

## Import VAST-TOOLS results tables
[http://matt.crg.eu/#get_vast](http://matt.crg.eu/#get_vast)
Check tables 
```bash
# check column names of data table
matt get_colnms DiffAS-Hsa14-hg19-dPSI15-range5_VCP.d0-vs-VCP.d7.tab

# print rows aligned
matt prnt_tab matt get_colnms DiffAS-Hsa14-hg19-dPSI15-range5_VCP.d0-vs-VCP.d7.tab
 -W 30 | less -s -

# check all different types of splicing events
matt col_uniq matt get_colnms DiffAS-Hsa14-hg19-dPSI15-range5_VCP.d0-vs-VCP.d7.tab
TYPE | matt prnt_tab - -w 9 | less -S -
```

Extract all intron retention events across samples and extract Gene IDs from GTD file used for alignment in VAST-TOOLS

`matt get_vast ~/working/oliver/projects/airals/splicing/raphaelle_vast_tools/vast_out/DiffAS-Hsa14-hg19-dPSI15-range5_VCP.d0-vs-VCP.d7.tab -minqab LOW -minqglob N -complex IR,IR-S,IR-C -a VCP.d7 -b VCP.d0 -gtf ~/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf -f gene_id > ir_events.tab`

## Define categories

Use this command def_cats for creating a new column with group IDs. Generally, rows of table corresponds to events, while columns correspond to their features. Hence, you define group IDs (one for each row) wrt. entries in the feature columns using the same constraint syntax as used with command [get_rows](http://matt.crg.eu/#get_rows).

This command is important as many analyses essentially compare some groups of items (e.g. groups of exons, introns, sequences, or more), hence, the user needs to define the groups to be compared.  
  
**Example:**  With  [t1.tab](http://matt.crg.eu/#TB)  being the table from above

`matt def_cats ir_events.tab GROUP2 'g1=START[0,5000000] STRAND]+[' 
     'g2=START[0,5000000] !STRAND]+[' 'g3=START[5000001,10000000]'`

will output a table with column GROUP2 with group IDs g1, g1, g2, g3, g1 categorising each micro exon in table t1.tab wrt. the defined constraints.

## Intron Features

Use get_ifeatures command to retrieve 50 features of interest for introns. Introns need to be described by a table with basic information like their genomic coordinates and a gene ID of genes the introns belong to. If the table does not yet contain gene IDs, you might use the command  [retr_geneids](http://matt.crg.eu/#retr_geneids)  for extracting gene IDs from any GTF file for given genomic events, like exons, introns, genes

![enter image description here](http://matt.crg.eu/graphics/ov_introns.png)

To get an overview of all the intron features
`matt get_ifeatures explain`

Write a table with columns containing the feature values & extracted sequence of intron, up/downstream exon & splice sites. If only one or a few features are of interest, users can apply  [get_cols](http://matt.crg.eu/#cmpr_exons)  and extract specific feature columns only.

For intron features help page:
`matt get_ifeatures help`

`matt get_ifeatures ir_events.tab START END SCAFFOLD STRAND GENEID ~/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf ~/working/oliver/genomes/sequences/human/GRCh38.primary_assembly.genome.fa Hsap -f gene_id > ifeatures.tab`



## Coverage for introns of interest
To perform the a focussed analysis of the 167 retained introns identified using VAST-tools, 

Run Sections A, B,  C "Import time effect" from `import_VASTOOLS.R` script located in `/home/camp/ziffo/working/oliver/scripts/intron_retention` 

Then run `get_relative_coverage_inteactive.R` script.  

Also note other R scripts:
`characterise_introns.R`

This approach uses a combination of threshold but accounts for the depth of coverage between different samples. Uses [MaxEntScan of 5' and 3' end](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) to identify splice sites. 

First run the global analysis to get the big picture of splicing events in VCP vs CTRL.  SVD & PCA analysis clustering by mutation. 
Then run detailed analysis according to specific genes or features. Multivariate analysis of the 2 VCP mutants. 

to be run in R which calculates a ratio of intron sequence coverage and surrounding exons. 

1. First check that the gene where the event is occurring is expressed. 
2. Then check the event exhibits a change over time of at least 10%. 
3. Finally inspect visually all selected IR events occurring at Day 7 NPC stage in VCP mutant and at Day 14 pMN stage in CTRL to end up with the list of 167 events to get a high confidence list.

Import the results obtained from VAST-tools (remember analysis in VCP & CTRL over time performed initially independently) 

To then get a value of IR across diverse data-sets I then wrote the custom code that computed the ratio between coverage of the intron versus average coverage of the neighbouring exons. Then select the events of interest. The `get_relative_coverage_interactive.R` script is located in `/home/camp/ziffo/working/oliver/scripts/intron_retention`





# DEXSeq

https://bioconductor.org/packages/release/bioc/html/DEXSeq.html
https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf
https://bioconductor.org/packages/release/bioc/manuals/DEXSeq/man/DEXSeq.pdf
https://genome.cshlp.org/content/22/10/2008

Measures differential exon usage (DEU) which indicates alternative splicing. DEU also measures alternative transcript start sites & polyadenylation sites (differential usage of exons at 5' and 3' boundary of transcripts). 

Calculates ratio = number of transcripts from the gene containing this exon / total number of transcripts from the gene.
Then compares this ratio between conditions assessing the strength of the fluctuations (dispersion). 

Steps:
1. Alignment to genome using splice aware aligner eg STAR
2. Counts using DEX Seq Python scripts (utilises HTSeq)

```bash
ml HTSeq
ml Python/2.7.15-GCCcore-7.3.0-bare
ml R

# use 2 Python scripts (dexseq_prepare_annotation.py & dexseq_count.py)
# load R environment with DEXSeq in
source activate rtest
R
library("DEXSeq")
#exit R
q()

# Create output folder
mkdir -p DEXSeq

### RUN ANNOTATION SCRIPT
# set GTF - Ensembl (gencode)
GTF=/home/camp/ziffo/working/oliver/genomes/annotation/gencode.v28.primary_assembly.annotation.gtf
# set GFF output
OUT=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_prepare_annotation.py

#run dexseq_prepare_annotation.py script
python $SCRIPT $GTF $OUT


### RUN COUNT SCRIPT
# set GFF input made in annotation script
GFF=/home/camp/ziffo/working/oliver/genomes/annotation/DEXSeq.homo_sapiens.GRCh38.gencode.v28.gff
SCRIPT=/camp/home/ziffo/R/x86_64-pc-linux-gnu-library/3.5/DEXSeq/python_scripts/dexseq_count.py
#set BAM input file
SAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/SRR5483796.sam
OUT=/home/camp/ziffo/working/oliver/projects/airals/splicing/DEXSeq/SRR5483796.txt

#run dexseq_count.py
python $SCRIPT $GFF $SAM $OUT

#to deactivate environment
source deactivate
```
Now switch to R to run the DEXSeq analysis.

```r
# make new Rproject & save it in the relevant CAMP Cluster splicing folder.
setwd("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/DEXSeq")
suppressPackageStartupMessages( library( "DEXSeq" ) )

### PREPARE DATA ###
# read in files
countFiles = list.files("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/DEXSeq/", pattern=".txt", full.names=TRUE) 
basename(countFiles)
flattenedFile = list.files("/Volumes/lab-luscomben/working/oliver/genomes/annotation/", pattern="gff", full.names=TRUE) 
basename(flattenedFile)

# create sample table: 1 row for each library. columns for file name & read counts, sample
sampleTable = data.frame(
  row.names = c( "SRR5483788", "SRR5483789", "SRR5483790", "SRR5483794", "SRR5483795", "SRR5483796" ), 
  condition = c("VCP", "VCP", "VCP", "CTRL", "CTRL", "CTRL" ) )
# check table
sampleTable
# create DEXSeqDataSet 
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
# see structure of dxd
sampleAnnotation( dxd )
colData(dxd)
# see first 5 rows from count data - 12 columns (first 6 = reads mapping to exons, last 6 = counts mapping to rest of exons from same gene)
head( counts(dxd), 5 )
#show details of exon bins annotation
head( rowRanges(dxd), 3 )

### DEXSeq ANALYSIS ###
#normalise for different sequencing depths between samples
dxd = estimateSizeFactors( dxd )
# Perform differential exon usage test in a single command (this by passes the steps below - takes 10 mins)
dxd = DEXSeq(dxd,
       fullModel=design(dxd),
       reducedModel = ~ sample + exon,
       BPPARAM=MulticoreParam(workers=1),
       fitExpToVar="condition")

# OPTIONAL: dispersion estimation - takes 30 mins!
dxd = estimateDispersions( dxd )
# plot the per-exon dispersion estimates vs mean normalised count
plotDispEsts( dxd )

# test for differential exon usage (DEU) between conditions
dxd = testForDEU( dxd )
# estimate relative exon usage fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

### RESULTS ###
# summarise the resuls
dxr1 = DEXSeqResults( dxd )
# description of each column in DEXSeqResults
mcols(dxr1)$description
# how many exons are significant (with false discovery rate 10% - can change this threshold to more stringent)
table ( dxr1$padj < 0.1 )
# how many genes are affected
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
# MA plot: log fold change vs avg normalised count per exon. Significant exons in red
plotMA( dxr1, cex=0.8 )
### Detailed overview of analysis results in HTML - see the geneIDs of differentially spliced genes
DEXSeqHTML( dxr1, FDR=0.1, color=c("#FF000080", "#0000FF80") )

### VISUALISATION ###
# visualise results for individual genes
plotDEXSeq( dxr1, "ENSG00000116560", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# visualise transcript models to see isoform regulation
plotDEXSeq( dxr1, "ENSG00000116560", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# visualise individual samples rather than model effect estimates.
plotDEXSeq( dxr1, "ENSG00000116560", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# remove overall changes from the plots
plotDEXSeq( dxr1, "ENSG00000116560", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
```

# rMATS
ml Python/2.7.15-GCCcore-7.3.0-bare
ml numpy
ml OpenBLAS
ml STAR
ml ScaLAPACK
ml GSL
ml GCC

[User Tutorial](http://rnaseq-mats.sourceforge.net/user_guide.htm)
http://rnaseq-mats.sourceforge.net/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4280593/

MATS is a computational tool to detect differential alternative splicing events from RNA-Seq data. The statistical model of MATS calculates the P-value and false discovery rate that the difference in the isoform ratio of a gene between 2 conditions exceeds a given user-defined threshold. From the RNA-Seq data, MATS can automatically detect and analyze alternative splicing events corresponding to all major types of alternative splicing patterns. MATS handles replicate RNA-Seq data from both paired and unpaired study design.

![enter image description here](http://rnaseq-mats.sourceforge.net/splicing.jpg)

# Perform GO analysis of gene list of differentially spliced genes

### Data Preparation

```bash
# convert splicing_diff.tab > .csv
perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < input.tab > output.csv
```
Open CSV in Excel
Copy the significant Gene Symbols into DAVID (see Gene Enrichment Chapter)

## topGO

### Prepare data from VAST tools splicing output

```r
library(tidyverse)
setwd("/Volumes/lab-luscomben/working/oliver/projects/airals/splicing/vast_tools/vast_out") #set working directory where splicing files are stored on cluster
diff_splicing <- read_tsv("diff.splicing_mv0.2.tab") #import differential splicing table .tab file
diff_splicing <- data.table(diff_splicing)

diff_splicing$pval <- as.numeric(diff_splicing$`E[dPsi]`)
diff_splicing$MV   <- as.numeric(diff_splicing$`MV[dPsi]_at_0.95`)

genelistUp <- factor( as.integer( diff_splicing$MV > .2 & diff_splicing$pval > 0) )
names(genelistUp) <- rownames(diff_splicing)
genelistDown <- factor( diff_splicing$MV > 0.2 & diff_splicing$pval < 0 )
names(genelistDown) <- rownames(diff_splicing)
```

### Run the topGO commands & make bar plots for BP, CC & MF (up & downregulated)

```r
### GO analysis
library(topGO)
library(GOstats)
library(goseq)
library(org.Mm.eg.db)
library(ggplot2)

### Test UPREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata      <- new( "topGOdata", ontology = "BP", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_BP_table   <- GenTable( myGOdata, goTestResults )
enrich.ur.bp  <- transform(UR_BP_table,result1 = as.numeric(result1))
enrich.ur.bp$value <- -log10(enrich.ur.bp$result1)
# Plot GO p-values as as bar plot
dat.ur.bp        <- enrich.ur.bp$value # the -log10(P-value)
names(dat.ur.bp) <- enrich.ur.bp$Term #the description of your GO term
barplot(height = dat.ur.bp, space = 0.5, horiz=T,las=1,font.size = 10)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.ur.cc <- transform(UR_CC_table,result1 = as.numeric(result1))
enrich.ur.cc$value <- -log10(enrich.ur.cc$result1)
# Plot GO p-values as as bar plot
dat.ur.cc        <- enrich.ur.cc$value # the -log10(P-value)
names(dat.ur.cc) <- enrich.ur.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.ur.cc,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistUp, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
UR_MF_table <- GenTable( myGOdata, goTestResults )
enrich.ur.mf <- transform(UR_MF_table,result1 = as.numeric(result1))
enrich.ur.mf$value <- -log10(enrich.ur.mf$result1)
### Plot GO p-values as as bar plot
dat.ur.mf        <- enrich.ur.mf$value # the -log10(P-value)
names(dat.ur.mf) <- enrich.ur.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.ur.mf,horiz=T,las=1, font.size = 20)

### TEST DOWNREGULATED GENES
#Test Biological Processes BP sub-ontology
myGOdata <- new( "topGOdata", ontology = "BP", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_BP_table <- GenTable( myGOdata, goTestResults )
enrich.dr.bp <- transform(DR_BP_table,result1 = as.numeric(result1))
enrich.dr.bp$value <- -log10(enrich.dr.bp$result1)
### Plot GO p-values as as bar plot
dat.dr.bp        <- enrich.dr.bp$value # the -log10(P-value)
names(dat.dr.bp) <- enrich.dr.bp$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.bp,horiz=T,las=1, font.size = 20)

#Test Cellular Compartment (CC) sub-ontology
myGOdata <- new( "topGOdata", ontology = "CC", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_CC_table <- GenTable( myGOdata, goTestResults )
enrich.dr.cc <- transform(DR_CC_table,result1 = as.numeric(result1))
enrich.dr.cc$value <- -log10(enrich.dr.cc$result1)
### Plot GO p-values as as bar plot
dat.dr.cc        <- enrich.dr.cc$value # the -log10(P-value)
names(dat.dr.cc) <- enrich.dr.cc$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.cc,horiz=T,las=1, font.size = 20)

#Test Molecular function (MF) sub-ontology
myGOdata <- new( "topGOdata", ontology = "MF", allGenes = genelistDown, nodeSize = 10,
                 annot = annFUN.org, mapping = "org.Hs.eg.db", ID="Entrez" )
goTestResults <- runTest( myGOdata, algorithm = "elim", statistic = "fisher" )
DR_MF_table <-GenTable( myGOdata, goTestResults )
enrich.dr.mf <- transform(DR_MF_table,result1 = as.numeric(result1))
enrich.dr.mf$value <- -log10(enrich.dr.mf$result1)
### Plot GO p-values as as bar plot
dat.dr.mf        <- enrich.dr.mf$value # the -log10(P-value)
names(dat.dr.mf) <- enrich.dr.mf$Term #the description of your GO term
par(mfrow=c(1,1),mar=c(3,20,3,3),cex=0.7)  # artificially set margins for barplot to follow (need large left hand margin for names)
barplot(height = dat.dr.mf,horiz=T,las=1, font.size = 20)
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE5MTg4MjU2OTUsNzQ5MTk2OTAzLDEwOD
MwNTczNTEsMjg1MDMxNTQ4LC0xNDI5NzA0NjA3LC0xODU3MzM2
ODEyLC0xOTIyOTYzMTMxLDgxNzg0NDIxNywtOTczODc4NDgyLC
04OTk1OTI2MDgsLTQwMjc5OTAxNiwtNzI0OTg0OTk5LDE3Njky
NDA3MTEsLTE5NDY4NzM4OTQsLTE3NTg1OTk5OSwxMTc4OTgxMD
IyLC0xNTMyNjY5MzYzLC0yMzUzNjY0MjAsLTU4NjgxMDE0Nywz
MjA0MDcxMjBdfQ==
-->