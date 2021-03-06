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
ml STAR

#set bam input
BAM=/home/camp/ziffo/working/oliver/projects/vcp_fractionation/alignment/roi_bam/*.bam
#set OUT directory
OUT=/camp/home/ziffo/working/oliver/projects/vcp_fractionation/alignment/roi_bam/BigWigs

for SAMPLE in $BAM
do
	sbatch -N 1 -c 4 --mem=24GB --wrap="STAR --runMode inputAlignmentsFromBAM --runThreadN 24 --inputBAMfile $SAMPLE --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix $SAMPLE"
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

Upload merged bigwig files to genome browser to download figures for publication.

1. Merge BAM files for each group
takes 10mins
```bash
ml SAMtools
ml STAR
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
sbatch -N 1 -c 4 --mem=24GB --wrap="STAR --runMode inputAlignmentsFromBAM --runThreadN 24 --inputBAMfile VCP_NPC_D7_merged.bam --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix VCP_NPC_D7_merged"

sbatch -N 1 -c 4 --mem=24GB --wrap="STAR --runMode inputAlignmentsFromBAM --runThreadN 24 --inputBAMfile VCP_iPSC_D0_merged.bam --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix VCP_iPSC_D0_merged"
```

3. Upload to [https://www.biodalliance.org/human38.html](https://www.biodalliance.org/human38.html)

Click + in top right of screen
Click Binary
Choose Files > select merged bigWig files
Edit the track - remove unwanted tracks, go to relevant genes

---
Alternatives are to upload to UCSC - couldnt get this to work as bigWigs too large.
[See page 13 here](https://www.france-bioinformatique.fr/sites/default/files/EBA/V3-2014/Nicolas_Servant_UCSC_TP_EBA-2014-10.pdf)
Login to UCSC browser: MyData > MySessions > Login [http://genome.ucsc.edu/cgi-bin/hgSession?hgsid=733417633_FrgAYw50hwa3PWBTCps2DLCXaSLK&hgS_doMainPage=1](http://genome.ucsc.edu/cgi-bin/hgSession?hgsid=733417633_FrgAYw50hwa3PWBTCps2DLCXaSLK&hgS_doMainPage=1)
Click Genome Browser in Menu.
Check reference genome is hg38
Enter gene of interest in search box > click track of interest
Click **Add custom track**, next to the bottom box click **upload**, on CAMP select the merged.bw file - allow time for upload

Or use **ggbio** in Rstudio
[http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf](http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf)

4. Download or screenshot

5. Finalise in Illustrator

**Sashimi Plots**
Visualise splice junctions & explore exon usage
![sashimi plot explained](http://miso.readthedocs.io/en/fastmiso/_images/sashimi-plot-example-annotated.png)
Bar graph height =  read coverage
Arcs = splice junctions
Numbers = number of reads that contain the respective splice junction.
IGV does not normalise for read number per sample in sashimi plots so dont overinterepret the read counts.

# UCSC genome browser

Follow instructions in this video: https://www.youtube.com/watch?v=UvHihNbyCh8

Create a public link to a UCSC genome browser session displaying the uploaded sequence tracks: https://genome.ucsc.edu/cgi-bin/hgGateway

Instructions are available here
http://genome.ucsc.edu/goldenPath/help/hgSessionHelp.html#NAR

## Create bigWig files

BAM > bigwigs

[https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
```bash
module load Anaconda3
python -m venv /camp/home/ziffo/home/envs 
source /camp/home/ziffo/home/envs/bin/activate  
python -m pip install deeptools

# For merged BAM files:
cd /camp/home/ziffo/home/projects/inter-neuron-bulk-rnaseq/alignment/STAR/merged
sbatch -N 1 -c 10 --mem=40G -t 48:00:00 --wrap="bamCoverage -b in_d18_ctrl.bam -o in_d18_ctrl.bw" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=bamCoverage --output=bamCoverage-%j.out --error=bamCoverage-%j.err
```

Convert GTF > BED > bigBed format




Move files to shareable location (Dropbox) and make shareable with others
```bash
cd /Users/ziffo/'Dropbox (The Francis Crick)'/Medical\ Files/Research/PhD/projects/inter-neuron-bulk-rnaseq/bigWigs
rsync -aP camp-ext:/camp/home/ziffo/home/projects/inter-neuron-bulk-rnaseq/alignment/STAR/merged/*.bw .
```
Create link > Copy link > paste into track lines script

## Create track lines settings script

```bash
track type=bigWig name="IN D18" description="interneuron day 18 control" color=0,0,255 visibility=2 bigDataUrl=https://www.dropbox.com/s/g7bond1pr88xqym/in_d18_ctrl.bw?dl=1
track type=bigWig name="IN D25" description="interneuron day 25 control" color=255,0,0 visibility=2 bigDataUrl=https://www.dropbox.com/s/jt18a30ik9qz5zk/in_d25_ctrl.bw?dl=1
track type=bigWig name="MN D18" description="motor neuron day 18 control" color=0,0,255 visibility=2 bigDataUrl=https://www.dropbox.com/s/jvyrath43s7kp33/mn_d18_ctrl.bw?dl=1
track type=bigWig name="MN D25" description="motor neuron day 25 control" color=255,0,0 visibility=2 bigDataUrl=https://www.dropbox.com/s/rsn2waae8m8mi6z/mn_d25_ctrl.bw?dl=1
track type=BED name="annotation" description="GENCODE annotation" color=0,255,0 visibility=3 bigDataUrl=https://www.dropbox.com/s/ke2ddj4zwl6bs74/Human.GRCh38.GENCODEv24.bed?dl=1
```
Ensure dropbox link is dl=1 at end

## Upload to UCSC
 
Login to UCSC

Reset user settings
http://genome.ucsc.edu/cgi-bin/cartReset

Configure genome browser
http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=968366871_acZUce1tkq5mADRGXsXv022MGtG8

Add custom tracks

In top box > paste the track lines settings into box > submit

Save as a new UCSC session:
http://genome.ucsc.edu/cgi-bin/hgSession?hgsid=968355787_90NOUDSzGcMJd8jgzNvA2TiOUEW2&hgS_doMainPage=1


After creating your custom tracks and viewing your data on the Genome Browser, you can save all of your tracks and settings to a snapshot of the Genome Browser called a  _session_. You can easily save a session by following these five steps:

1.  **Configure the Genome Browser to your preference**  
    Make sure the display of your custom tracks is to your liking on the Genome Browser.
2.  **Navigate to the sessions page**  
    Once you are satisfied with the display, go to the  [My Sessions](http://genome.ucsc.edu/cgi-bin/hgSession)  page by either:
    -   Going to  _**My Data**_  ->  _**My Sessions**_  from the navigation bar.
    -   Using the "_**s**  then  **s**_" keyboard shortcut when viewing the main page of the  [Genome Browser](http://genome.ucsc.edu/cgi-bin/hgTracks).
3.  **Login to the UCSC Genome Browser**  
    You must sign in to be able to save named sessions which will then be displayed with Browser and Email links.
4.  **Save your session**  
    Go to the  **Save Settings**  section and in the  _Save current settings as named session_  text box, enter a name for your session. When saving the session, be sure to have the "_allow this session to be loaded by others_" option checked and then click  _**submit**_.
5.  **Edit the session description**  
    Once the session is created, you can click the  _**details**_  button to add a description to the session. If you eventually make a Public Session, and provide a detailed description, anyone can find it by searching for terms you share. For example, navigate to the  [Public Sessions](http://genome.ucsc.edu/cgi-bin/hgPublicSessions)  page and search "NAR" to see some example sessions.

<!--stackedit_data:
eyJoaXN0b3J5IjpbLTQ4NDI3NDExMyw1ODgxMzI0MTMsNjQ5Mj
k5NCw0MTU0NDE5NjUsLTE5NDEzODM4NTUsMjQ1ODE2NzgxLDg4
MTQwODQwMSwtMTkxNjEwNjQxMiw3NTIzOTI2OTcsLTE5MjUxMD
Q4NzAsLTI1Nzk1OTY2OCwxMTAyNjAyNDc2LDYwMjIyMzgyMywx
MDYxMzY3MjQ3LDE2ODk2MzY1NjcsNTA2OTU3MDIyLDYwODk0Mz
kzMCwxMzcwMDcxODYsMTM5MzI1NTY5MywxNjc5MjE1MjY4XX0=

-->