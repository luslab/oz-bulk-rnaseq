


> # Alternative Expression Analysis
By infer structural infromation about the transcript 
Infer the strand by examining splice site spanning reads

# Tools
Cuffdiff [14], DEXSeq [192], ALEXA-seq [3], IUTA [193], FineSplice [194], PennSeq [195], FlipFlop [196], SNPlice [197], spliceR [198], GESS [102], RNASeq-MATS [199], SplicingCompass [200], DiffSplice [201], SigFuge [202], SUPPA [bioRXiv], CLASS [bioRXiv], SplAdder [bioRXiv], SplicePie [203].

# Gene Isoform counting

Programmes to quantify isoforms:
`Cufflinks`
`RSEM`
`eXpress`

Gene isoforms are mRNA produced from the same locus but with diferent protein codeing sequences:

5 modes of alternative splicing are recognized:

1.  Exon skipping or cassette exon.
2.  Mutually exclusive exons.
3.  Alternative donor (5') site.
4.  Alternative acceptor (3') site.
5.  Intron retention.

![Alternative Splicing](https://en.wikipedia.org/wiki/Protein_isoform#/media/File:Alternative_splicing.jpg)

**Preparing an annotation:**
To **assess differential expression of exons**, create an annotation file where overlapping exons of different isoforms are split before running featureCounts. Use `dexseq_prepare_annotation.py` script of DEXSeq package or `QoRTs`.

@Raphaelle used: 
 **Splicing analysis**  
   
-   I first used  [VAST-tools](https://github.com/vastgroup/vast-tools) which performs alignment for you. So basically you submit your fastq files directly. Have a look at the GitHub vignette as it is rather complete however please do not hesitate to contact me if you want help with shell scripting as you will need to run this as a loop. Or Nobby will certainly be happy to help on CAMP (I am working from UCL cluster and have never logged onto CAMP).
-   Then to perform the more focussed analysis on the 167 retained introns, which I identified using VASt-tools, I wrote a script in R which basically obtain the coverage for intronic sequences of interest and surrounding exons and then compute the ratio. As input I use the BAM files.

  **3' UTR isoforms analysis**
-   I have written an entire pipeline for this which I can explain and share with you scripts when needed. But basically I extract genome-wide coverage using bedtools, then extract regions of continuous coverage along genome, then intersect these with Ensembl annotated regions, extend 3' UTR. Finally to annotate all alternative 3' UTR isoforms I then run an algo which identifies shifts in coverage along 3' UTR which are expected to occur at PAS sites.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEzMzU0NzUwMCwtNTQyMzA4MzY5XX0=
-->