---


---

<blockquote>
<h1 id="quality-control">Quality Control</h1>
</blockquote>
<p><a href="https://osca.bioconductor.org/quality-control.html">https://osca.bioconductor.org/quality-control.html</a><br>
Low-quality libraries in scRNA-seq data can arise from:</p>
<ul>
<li>cell damage during dissociation</li>
<li>failure in library preparation (e.g., inefficient reverse transcription or PCR amplification).</li>
</ul>
<p>These usually manifest as “cells” with:</p>
<ul>
<li>low total counts</li>
<li>few expressed genes</li>
<li>high mitochondrial or spike-in proportions.</li>
</ul>
<p>Low-quality libraries are problematic as they mislead results:</p>
<ul>
<li>
<p>They form their own distinct <strong>cluster</strong>(s), complicating interpretation of results. Low-quality libraries generated from different cell types can cluster together based on similarities in the damage-induced expression profiles, creating <strong>artificial intermediate states</strong> or trajectories between otherwise distinct subpopulations. Very small libraries can form their own clusters due to shifts in the mean upon transformation</p>
</li>
<li>
<p>They distort the characterization of population heterogeneity during variance estimation or principal components analysis. The first few <strong>principal components will capture differences in quality rather than biology</strong>, reducing effectiveness of dimensionality reduction. Genes with the largest variances will be driven by differences between low- and high-quality cells.</p>
</li>
<li>
<p>Low-quality libraries with very low counts where scaling normalization <strong>inflates the apparent variance of genes</strong> that happen to have a non-zero count in those libraries. Genes that appear to be strongly “upregulated” are due to <strong>aggressive scaling</strong> to normalize for small library sizes. Contaminating transcripts (e.g., from the ambient solution) that are present in all libraries at low but constant levels. Increased scaling in low-quality libraries transforms small counts for these transcripts in large normalized expression values, resulting in apparent upregulation compared to other cells. This can be misleading as the affected genes are often biologically sensible but are actually expressed in another subpopulation.</p>
</li>
</ul>
<p>To mitigate this, we need to <strong>remove these cells at the start</strong> of the analysis. This step is commonly referred to as quality control (QC) on the cells. We use “library” and “cell” interchangeably, though the distinction will become important when dealing with droplet-based data.</p>
<h1 id="qc-metrics">QC metrics</h1>
<p><strong>Library size</strong> = total sum of counts across all relevant features ( = endogenous genes) for each cell. Cells with small library sizes are of low quality as the RNA has been lost at some point during library preparation (either cell lysis or inefficient cDNA capture and amplification).</p>
<p><strong>Number of expressed features in each cell</strong> = number of endogenous genes with non-zero counts for that cell. Any cell with very few expressed genes is likely to be of poor quality.</p>
<p><strong>Proportion of reads mapped to spike-in transcripts</strong> is calculated relative to the total count across all features (including spike-ins) for each cell. As the same amount of spike-in RNA should have been added to each cell, any enrichment in spike-in counts is symptomatic of loss of endogenous RNA. High proportions are indicative of poor-quality cells where endogenous RNA has been lost.</p>
<p><strong>Proportion of reads mapped to genes in the mitochondrial genome</strong> can be used when spike-ins are not available. High proportions = poor-quality cells because of <strong>loss of cytoplasmic RNA from perforated cells</strong>. In modest damage the holes in the cell membrane permit efflux of individual transcript molecules but are too small to allow mitochondria to escape, leading to a relative enrichment of mitochondrial transcripts.</p>
<p>Outliers can also be identified from the <strong>gene expression profiles,</strong> rather than QC metrics. This is a risky strategy as it can remove high-quality cells in rare populations.</p>
<h2 id="scatter-package">Scatter package</h2>
<p>For each cell, we calculate these QC metrics using the  <code>perCellQCMetrics()</code>  function from the  <em><a href="https://bioconductor.org/packages/3.10/scater">scater</a></em>  package.</p>
<p><code>sum</code>  column contains the total count for each cell<br>
<code>detected</code>  column contains the number of detected genes<br>
<code>subsets_Mito_percent</code>  contains the percentage of reads mapped to mitochondrial transcripts (based on Ensembl annotation)<br>
<code>altexps_ERCC_percent</code>  contains the percentage of reads mapped to ERCC transcripts.</p>
<h1 id="identify-low-quality-cells">Identify low quality cells</h1>
<h2 id="fixed-thresholds">Fixed thresholds</h2>
<p>The simplest approach to identifying low-quality cells is to apply thresholds on the QC metrics. e.g. consider cells to be low quality if they have:</p>
<ul>
<li>library sizes below 100,000 reads</li>
<li>express fewer than 5,000 genes</li>
<li>spike-in proportions above 10%</li>
<li>have mitochondrial proportions above 10%</li>
</ul>
<p>While simple, this strategy <strong>requires considerable experience</strong> to determine appropriate thresholds. Thresholds for read count-based data are simply not applicable for UMI-based data, and vice versa.</p>
<h2 id="adaptive-thresholds">Adaptive thresholds</h2>
<p>Outlier detection assumes that most cells consists of high-quality cells. Then identify cells that are outliers for the various QC metrics, based on the median absolute deviation (MAD) from the median value of each metric across all cells.</p>
<p>An outlier = &gt; 3 MADs from the median in the “problematic” direction. This filter will retain ~99% of non-outlier values that follow a normal distribution.</p>
<p>A log-transformation is used to improve resolution at small values and ensures the threshold is not a negative value, which would be meaningless. These metrics can exhibit a heavy right tail, and the log-transformation makes the distribution seem more normal to justify the 99% rationale.</p>
<p>If most cells are of (unacceptably) low quality, the adaptive thresholds will fail as they cannot remove the majority of cells.</p>
<p>The assumption that QC metrics are independent of the biological state of each cell is violated in highly heterogeneous cell populations where cell types that naturally have less RNA or more mitochondria are more likely to be considered outliers and removed, even if they are of high quality. MAD mitigates this to some extent by accounting for biological variability in the QC metrics. A heterogeneous population has higher variability in the metrics among high-quality cells, increasing the MAD and reducing the chance of incorrectly removing particular cell types (at the cost of reducing power to remove low-quality cells).</p>
<h2 id="technical-replicates">Technical replicates</h2>
<p>Technical repicates with different sequencing depth should apply adaptive thresholds to each batch separately. Computing medians and MADs on the mixture distribution with samples from multiple batches will alter the mean &amp; MADs. If sequencing coverage is lower in one batch, it will drag down the median and inflate the MAD. This reduces the suitability of the adaptive threshold for the other batches.</p>
<p>If each batch is represented by its own  <code>SingleCellExperiment</code>,  <code>isOutlier()</code>  can be directly applied to each batch separately. However, if cells from all batches have been merged into a single  <code>SingleCellExperiment</code>, the  <code>batch=</code>  argument should be used to ensure that outliers are identified  <em>within</em>  each batch. This allows  <code>isOutlier()</code>  to accommodate systematic differences in the QC metrics across batches.</p>
<p><code>batch=</code> makes the assumption that most cells in each batch are of high quality. If an entire batch failed, outlier detection will not be able to act as an appropriate QC filter for that batch. This inflates the median and MAD within those batches, resulting in a failure to remove the assumed low-quality cells. In such cases, it is better to either not use <code>batch=</code> or to apply a custom filter to the problematic batches.</p>
<h1 id="check-diagnostic-plots">Check Diagnostic Plots</h1>
<ol>
<li>Total counts</li>
<li>Detected features (genes)</li>
<li>Mitochondrial %</li>
</ol>
<p>Also plot the proportion of mitochondrial counts against some of the other QC metrics. The aim is to confirm that there are no cells with both large total counts and large mitochondrial counts, to ensure that we are not inadvertently removing high-quality cells that happen to be highly metabolically active</p>
<p>It is crucial to inspect the distributions of QC metrics. Ideally we would see normal distributions that justify the 3 MAD threshold used in outlier detection.</p>
<p>A large proportion of cells in another mode suggests that the QC metrics might be correlated with some biological state, potentially leading to the loss of distinct cell types during filtering.</p>
<p>Batches with systematically poor values for any metric can also be quickly identified for further troubleshooting or outright removal.</p>
<h1 id="identify-empty-droplets">Identify empty droplets</h1>
<p>With droplet-based data we dont know whether a particular library (i.e., cell barcode) corresponds to cell-containing or empty droplets. We need to call cells from empty droplets based on the <strong>observed expression profiles</strong>.</p>
<p>This is not straightforward as empty droplets can contain ambient (i.e., extracellular) RNA that can be captured and sequenced, resulting in non-zero counts for libraries that do not contain any cell.</p>
<p>The distribution of total counts exhibits a sharp transition between barcodes with large and small total counts which probably corresponds to cell-containing and empty droplets. A simple approach is to apply a threshold on the total count to only retain those barcodes with large totals. However, this unnecessarily discards libraries derived from cell types with low RNA content.</p>
<p>Use <code>emptyDrops()</code> to test whether the <strong>expression profile for each cell barcode is significantly different from the ambient RNA pool.</strong> Any significant deviation indicates the barcode corresponds to a cell-containing droplet. This allows discrimination between empty droplets and cells with little RNA, both of which would have similar total counts in. Call cells at a false discovery rate (FDR) of 0.1%, i.e. no more than 0.1% of our called barcodes should be empty droplets on average.</p>
<p><code>emptyDrops()</code>  uses Monte Carlo simulations to compute  p-values for the multinomial sampling transcripts from the ambient pool. The number of Monte Carlo iterations determines the lower bound for the  p-values. The  <code>Limited</code>  field in the output indicates whether or not the computed  p-value for a particular barcode is bounded by the number of iterations. If any non-significant barcodes are  <code>TRUE</code>  for  <code>Limited</code>, we may need to increase the number of iterations. A larger number of iterations will result in a lower p-value for these barcodes, which may allow them to be detected after correcting for multiple testing.</p>
<p><code>emptyDrops()</code> assumes that barcodes with low total UMI counts are empty droplets. The null hypothesis should be true for all of these barcodes. We can check whether the hypothesis testing procedure holds its size by examining the distribution of p-values for low-total barcodes with <code>test.ambient=TRUE</code>. Ideally, the distribution should be close to uniform. Large peaks near zero indicate that barcodes with total counts below <code>lower</code> are not all ambient in origin. This can be resolved by decreasing <code>lower</code> further to ensure that barcodes corresponding to droplets with very small cells are not used to estimate the ambient profile.</p>
<p>Once satisfied with <code>emptyDrops()</code>, subset the <code>SingleCellExperiment</code> object to retain only detected cells. Discerning readers will notice the use of <code>which()</code>, which conveniently removes the <code>NA</code>s prior to the subsetting.</p>
<p>We do not attempt to remove the ambient contamination from each library. Accurate quantification of the contamination rate in each cell is difficult as it generally requires some prior biological knowledge about genes that are expected to have mutually exclusive expression profiles <em>and</em> are highly abundant in the ambient solution. Fortunately, ambient contamination usually has little effect on the downstream conclusions for routine analyses; cell type identities are usually easy enough to determine from the affected genes, notwithstanding a (mostly harmless) low background level of expression for marker genes that should be unique to a cell type.</p>
<p>While <code>emptyDrops()</code> distinguishes cells from empty droplets, it makes no statement about the quality of the cells. It is possible for droplets to contain damaged or dying cells, which need to be removed prior to downstream analysis. Filtering on the mitochondrial proportion provides the most additional benefit in this situation, provided that we check that we are not removing a subpopulation of metabolically active cells.</p>
<p><code>emptyDrops()</code>  already removes cells with very low library sizes or (by association) low numbers of expressed genes. Thus, further filtering on these metrics is not strictly necessary. It may still be desirable to filter on both of these metrics to remove non-empty droplets containing cell fragments or stripped nuclei that were not caught by the mitochondrial filter. However, this should be weighed against the risk of losing genuine cell types.</p>
<p><em>CellRanger</em>  version 3 automatically performs cell calling using an algorithm similar to  <code>emptyDrops()</code>. If we had started our analysis with the  <strong>filtered</strong>  count matrix, we could go straight to computing other QC metrics. We do not need to run  <code>emptyDrops()</code> and indeed, attempting to do so would lead to nonsensical results if not outright software errors.</p>
<h1 id="remove-low-quality-cells">Remove low quality cells</h1>
<p>Once low-quality cells have been identified, we can choose to either remove them or mark them. Removal is the most straightforward option and is achieved by subsetting the  <code>SingleCellExperiment</code>  by column.</p>
<p>The biggest concern during QC is whether an entire cell type is inadvertently discarded. There is always some risk of this occurring as the QC metrics are never fully independent of biological state.</p>
<p>We can diagnose cell type loss by looking for systematic differences in gene expression between the discarded and retained cells. To demonstrate, we compute the average count across the discarded and retained pools in the data set, and we compute the log-fold change between the pool averages. If the discarded pool is enriched for a certain cell type, we should observe increased expression of the corresponding marker genes.</p>
<p>If we suspect that cell types have been incorrectly discarded, the most direct solution is to relax the QC filters for metrics that are associated with genuine biological differences. Of course, this increases the risk of retaining more low-quality cells. The logical endpoint of this line of reasoning is to avoid filtering altogether.</p>
<p>The true technical quality of a cell may also be correlated with its type. (This differs from a correlation between the cell type and the QC metrics, as the latter are our imperfect proxies for quality.) This can arise if some cell types are not amenable to dissociation or microfluidics handling during the scRNA-seq protocol. In such cases, it is possible to “correctly” discard an entire cell type during QC if all of its cells are damaged. Indeed, concerns over the computational removal of cell types during QC are probably minor compared to losses in the experimental protocol.</p>
<h1 id="mark-low-quality-cells">Mark low quality cells</h1>
<p>The other option is to simply mark the low-quality cells as such and retain them in the downstream analysis. The aim here is to allow clusters of low-quality cells to form, then to identify and ignore such clusters during interpretation of the results. This approach avoids discarding cell types that have poor values for the QC metrics. This gives an opportunity to decide whether a cluster of such cells represents a genuine biological state.</p>
<p>The downside is that it shifts the burden of QC to the interpretation of the clusters, which is already the bottleneck in scRNA-seq data analysis (Chapters  <a href="https://osca.bioconductor.org/clustering.html#clustering">10</a>,  <a href="https://osca.bioconductor.org/marker-detection.html#marker-detection">11</a>  and  <a href="https://osca.bioconductor.org/cell-type-annotation.html#cell-type-annotation">12</a>). Indeed, if we do not trust the QC metrics, we would have to distinguish between genuine cell types and low-quality cells based only on marker genes, and this is not always easy due to the tendency of the latter to “express” interesting genes (Section  <a href="https://osca.bioconductor.org/quality-control.html#quality-control-motivation">6.1</a>).</p>
<p>Retention of low-quality cells also compromises the accuracy of the variance modelling, requiring, e.g., use of more PCs to offset the fact that the early PCs are driven by differences between low-quality and other cells.</p>
<p>For routine analyses, perform removal by default to avoid complications from low-quality cells. This allows most of the population structure to be characterized with no - or, at least, fewer - concerns about its validity. Once the initial analysis is done, and if there are any concerns about discarded cell types, a more thorough re-analysis can be performed where the low-quality cells are only marked. This recovers cell types with low RNA content, high mitochondrial proportions, etc. that only need to be interpreted insofar as they “fill the gaps” in the initial analysis.</p>

