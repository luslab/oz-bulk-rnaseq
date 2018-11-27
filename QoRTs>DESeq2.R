library(QoRTs)

res <- read.qc.results.data("/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC/", 
                            decoder.files = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.byUID.txt", 
                            calc.DESeq2 = TRUE, calc.edgeR = TRUE); 

# Once you have read in the QC data, you can build all sorts of plots.
# EXAMPLE 1: The makeMultiPlot.all can be used to automatically generate a full battery of multi-plot figures: 
makeMultiPlot.all(res, 
                  outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPlots/", 
                  plot.device.name = "png"); 

#EXAMPLE 2: Some users may find the large png files difficult to read. QoRTs offers multi-page pdf reports as an alternative, simply by using the plot.device.name parameter: 
makeMultiPlot.all(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/summaryPDFs/",
                  plot.device.name = "pdf"); 

#EXAMPLE 3: To print all the basic plots as seperate pngs, use the command: 
makeMultiPlot.basic(res, outfile.dir = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/basicPlots/", 
                    separatePlots = TRUE);

# Extract size factors. QoRTs generates these to normalise all samples to a comparable scale allowing downstream comparison with DESeq2 or edgeR
get.size.factors(res, outfile = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/sizeFactors.GEO.txt");


suppressPackageStartupMessages(library(DESeq2))

decoder.bySample <- read.table("/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/QoRTs/decoder.bySample.txt", 
                               header=T,stringsAsFactors=F); 

directory <- "/Volumes/lab-luscomben/working/oliver/projects/airals/alignment/D7_samples/trimmed_filtered_depleted/alignment_QC"; 
sampleFiles <- paste0(decoder.bySample$qc.data.dir, "/QC.geneCounts.formatted.for.DESeq.txt.gz" ); 

sampleCondition <- decoder.bySample$group.ID; 
sampleName <- decoder.bySample$sample.ID; 
sampleTable <- data.frame(sampleName = sampleName, 
                          fileName = sampleFiles, 
                          condition = sampleCondition); 

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                  directory = directory, 
                                  design = ~ condition); 
dds

dds <- DESeq(dds);
res <- results(dds);
res;

# Sort the results data frame by the padj and foldChange columns. 
sorted = res[with(res, order(padj,  -log2FoldChange)),  ]  
# Turn it into a dataframe to have proper column names. 
sorted.df = data.frame("id"=rownames(sorted),sorted)
#write the table out:
write.table(sorted.df, file = "/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2/DESeq2.results.txt", sep="\t", col.names=NA, quote=FALSE);

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="/Volumes/lab-luscomben/working/oliver/projects/airals/expression/D7_samples/DESeq2/norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)