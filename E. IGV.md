># Genome Browser
- Check results visually to ensure reads align to expected regions without excess mismatches.
- Genome browsers give a linear track representing the forward stand of the genome. left = 5'. right = 3'
- Visualise in line orientated formats (fasta, bed, gff, SAM/BAM)
- genomic features are drawn over the linear track in **glyphs** (pictogram)
	- Horizontal intervals = directions, genes, alignments
	- Values over intervals = coverages, probabilities
	- Attributes at locations = mutations, deletions, junctions
- To create publication genome browser shots, merge bedgraphs for each group, then upload tracks to [UCSC genome browser](http://genome.ucsc.edu/index.html) > download

### `samtools tview`
- simplest genome browser. can visualise any DAM file `samtools tview --reference reference_genome.fa filename.bam`

## Standalone Genome Browsers
- [Integrative Genomics Viewer IGV](http://software.broadinstitute.org/software/igv/book/export/html/6)  by the Broad Institute
- [Integrated Genome Browser IGB](https://bioviz.org/) by University of North Carolina
- [Artemis](https://www.sanger.ac.uk/science/tools/artemis) by the Wellcome Sanger Institute
- [Tablet](https://ics.hutton.ac.uk/tablet/) by James Hutton Institute
- [SeqMonk](https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/) for high throughput RNA-seq
- [iobio](http://iobio.io/) for real-time genomic analysis

#### Online Genome Browsers
- [Ensembl](http://useast.ensembl.org/index.html)
- [UCSC](https://genome.ucsc.edu/)

# Convert BAM > Big Wig file coverage tracks
Convert BAM > BigWig file using BED tools & BED graph. Can then import the BigWig file into IGV & UCSC genome browsers. 
STAR creates a wiggle track (raw bigwig file). 
They only show the coverage and not the individual reads. They are binary formatted and so a much smaller files - prevents IGV crashing.

https://software.broadinstitute.org/software/igv/bigwig
https://github.com/YangLab/bamTobw

Use bamCoverage to convert BAM > BW files
```bash
ml SAMtools
ml IGVTools

#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/airals/alignment/D0_samples/*Aligned.sortedByCoord.out.bam
#set OUT directory
OUT=/home/camp/ziffo/working/oliver/projects/airals/alignment/D0_samples/

for SAMPLE in $BAM
do
	SRRID=`echo $SAMPLE | grep -E -o 'SRR[0-9]+'`
	sbatch -N 1 -c 4 --mem=24GB --wrap="bamCoverage -b $SAMPLE -o $OUT_$SRRID.bw"
done
```

# IGV
Best resources are the [IVG mannual](http://software.broadinstitute.org/software/igv/userguide) and [youtube videos](https://www.youtube.com/results?search_query=integrative+genome+viewer)
https://www.biostarhandbook.com/visualize/igv.html

1. Run IGV on local computer and mount CAMP. [Set Java 8 as default](https://stackoverflow.com/questions/46513639/how-to-downgrade-java-from-9-to-8-on-a-macos-eclipse-is-not-running-with-java-9) since IGV doesnt work with Java 10
```bash
# On local terminal move to IGV in bin
cd ~/bin/IGV_2.4.14/lib
# run IGV via command line on local terminal
java -Xmx750m -jar igv.jar
```
2. Set reference genome to Human (hg38) top left box.
3. 
4. Click File load from file > click Desktop > mount CAMP locally > click relevant BigWig (.bw) or BAM files (to see individual read tracks)  files (load multiple at once).

To visualise on IGV its easier to generate TDF files which are much lighter. This is useful if want to add more data-sets later. To generate TDF files first generate Bedgraph coverage files, then sort and then create the tdf file. Create 3 different coverage files: positive, negative strands and total. As the data is stranded it is better to look at both strands separately. Run the code in file: PE_strandedBedGraph.sh

4. Rename BAM files: right click file name in left hand column.

5. Go to the Genomic Location of interest: search gene name or gene coordinates e.g. SFPQ or chr1:35,176,378-35,193,158 in top middle box > Go. Find Genomic Locations using google search eg https://www.genecards.org/cgi-bin/carddisp.pl?gene=SFPQ
Mark the region of interest: Regions > Region Navigator > Add

6. Zoom in & out using top right zoomer or +/-. IGV will ask you to zoom in to see alignment data. This is because it has a default of 30 kb memory. You can change this to a higher value to see alignment data from a further out zoom but this will slow down the processing speed when scrolling across the genome. If you have deep coverage files then keep memory low to avoid it slowing down processing.

For each BAM file there are 2 tracks. Top track shows summary distribution of the coverage of the exonic islands separated by spice junctions (introns). Bottom tracks show all the individual sequence read alignments piled up.

7. Compare exon usage in the top tracks to transcript genome (bottom row): make the rows for each of the BAM files smaller and right click on the reference genome at the bottom row (in blue) > expand to show all the previously annotated differential splicing.

Bases that dont match the reference sequences are highlighted by colour. The deeper the shade of grey the more confidence you can have that the sequence was aligned correctly. White means no confidence alignment.

![IGV RNA-seq specific view](https://lh3.googleusercontent.com/h7PbqBtb3kHxxevIpjvKJUAd451K0UFOoACMogIZzUhVVMz-_AqRnjSYsNpmhYeCbct9ikfaZU8-Yg "IGV RNA-seq specific view")

IGV is used only to validate & confirm analysis results.  Use it to explore large genomic datasets. It is not good for the primary analysis.
![IGV of SFPQ D7 NPC samples](https://lh3.googleusercontent.com/r8Ph08oRuLWUBmnc6gbEyX5Rg3iBEkGhNmmNTHqTr7J01dtwdBGIdAqYJ2BMNlLcVIyYxPbn0QEhTQ "IGV of SFPQ D7 NPC samples")

Top row = chromosome. red bar is location. blue lines mid-section refer to transcripts binding with more = higher peak. bottom section = reference genomes.
Coverage line = quick identification of highly covered regions
Grey boxes = aligned reads. 
Gaps between reads with horizontal grey line = introns
Blue lines within reads = insertions
Red lines within reads = deletions
Blue boxes = reference genome

# Create Publication Genome Browser view
[Short introduction to UCSC genome browser PDF](https://www.france-bioinformatique.fr/sites/default/files/EBA/V3-2014/Nicolas_Servant_UCSC_TP_EBA-2014-10.pdf)
[http://genome.ucsc.edu/training/index.html](http://genome.ucsc.edu/training/index.html)
[https://www.youtube.com/channel/UCQnUJepyNOw0p8s2otX4RYQ](https://www.youtube.com/channel/UCQnUJepyNOw0p8s2otX4RYQ)
[https://www.sciencedirect.com/science/article/pii/S0888754308000451?via%3Dihub](https://www.sciencedirect.com/science/article/pii/S0888754308000451?via%3Dihub)
[http://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html](http://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html)

Upload merged bigwig files to UCSC genome browser to download figures for publication.

1. Merge BAM files for each group
takes 10mins
```bash
ml SAMtools
# VCP NPC D7
sbatch -N 1 -c 4 --mem=24GB --wrap="samtools merge VCP_NPC_D7_merged.bam SRR5483788_Aligned.sortedByCoord.out.bam SRR5483789_Aligned.sortedByCoord.out.bam SRR5483790_Aligned.sortedByCoord.out.bam"
# VCP iPSC D0
sbatch -N 1 -c 4 --mem=24GB --wrap="samtools merge VCP_iPSC_D0_merged.bam SRR5483800Aligned.sortedByCoord.out.bam SRR5483801Aligned.sortedByCoord.out.bam SRR5483802Aligned.sortedByCoord.out.bam"

# index merged BAMs
sbatch -N 1 -c 8 --mem 40 --wrap="samtools index VCP_NPC_D7_merged.bam"
sbatch -N 1 -c 8 --mem 40 --wrap="samtools index VCP_iPSC_D0_merged.bam"
```
2. Convert merged BAMs > BigWig & normalise for read depth
```bash
sbatch -N 1 -c 4 --mem=24GB --wrap="bamCoverage -b VCP_NPC_D7_merged.bam -o VCP_NPC_D7_merged.bw --normalizeUsing RPKM"
sbatch -N 1 -c 4 --mem=24GB --wrap="bamCoverage -b VCP_iPSC_D0_merged.bam -o VCP_iPSC_D0_merged.bw --normalizeUsing RPKM"
```
3. Upload to UCSC genome browser
[See page 13 here](https://www.france-bioinformatique.fr/sites/default/files/EBA/V3-2014/Nicolas_Servant_UCSC_TP_EBA-2014-10.pdf)

Login to UCSC browser: MyData > MySessions > Login [http://genome.ucsc.edu/cgi-bin/hgSession?hgsid=733417633_FrgAYw50hwa3PWBTCps2DLCXaSLK&hgS_doMainPage=1](http://genome.ucsc.edu/cgi-bin/hgSession?hgsid=733417633_FrgAYw50hwa3PWBTCps2DLCXaSLK&hgS_doMainPage=1)

Click Genome Browser in Menu.
Check reference genome is hg38
Enter gene of interest in search box > click track of interest

Click **Add custom track**, next to the bottom box click **upload**, on CAMP select the merged.bw file - allow time for upload



4. Optimise browser appearance
Change views of each track to "pack, denser or hide"

5. Download
6. Finalise in Illustrator

**Sashimi Plots**
Visualise splice junctions & explore exon usage
![sashimi plot explained](http://miso.readthedocs.io/en/fastmiso/_images/sashimi-plot-example-annotated.png)
Bar graph height =  read coverage
Arcs = splice junctions
Numbers = number of reads that contain the respective splice junction.
IGV does not normalise for read number per sample in sashimi plots so dont overinterepret the read counts.
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTA2MTM2NzI0NywxNjg5NjM2NTY3LDUwNj
k1NzAyMiw2MDg5NDM5MzAsMTM3MDA3MTg2LDEzOTMyNTU2OTMs
MTY3OTIxNTI2OCwxMjYwMTgxMTg0LC0xMjg4NTYxMTk1LC0xOD
M0MDI1NDI0LC0xMDQ0Njg1ODY1LDEzNTM5MTc4MjMsLTM4OTE0
MDI2OCw3NTc1MzExMDIsOTE4OTE5NTUyLDIxMTkxMDQ1NSwxOD
A5MzY5MTc3LDI4MzA4NzA0LDU3NDM0NzA4OSwtMTg1NzExOTU1
N119
-->