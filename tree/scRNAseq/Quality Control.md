---


---

<blockquote>
<h1 id="quality-control">Quality Control</h1>
</blockquote>
<p>Low-quality libraries in scRNA-seq data can arise from:</p>
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
<p>These low-quality libraries are problematic as they contribute to misleading results in downstream analyses:</p>
<ul>
<li>
<p>They form their own distinct <strong>cluster</strong>(s), complicating interpretation of results. Low-quality libraries generated from different cell types can cluster together based on similarities in the damage-induced expression profiles, creating <strong>artificial intermediate states</strong> or trajectories between otherwise distinct subpopulations. Very small libraries can form their own clusters due to shifts in the mean upon transformation</p>
</li>
<li>
<p>They distort the characterization of population heterogeneity during variance estimation or principal components analysis. The first few principal components will capture differences in quality rather than biology, reducing the effectiveness of dimensionality reduction. Similarly, genes with the largest variances will be driven by differences between low- and high-quality cells. The most obvious example involves low-quality libraries with very low counts where scaling normalization inflates the apparent variance of genes that happen to have a non-zero count in those libraries.</p>
</li>
<li>
<p>They contain genes that appear to be strongly “upregulated” due to aggressive scaling to normalize for small library sizes. This is most problematic for contaminating transcripts (e.g., from the ambient solution) that are present in all libraries at low but constant levels. Increased scaling in low-quality libraries transforms small counts for these transcripts in large normalized expression values, resulting in apparent upregulation compared to other cells. This can be misleading as the affected genes are often biologically sensible but are actually expressed in another subpopulation.</p>
</li>
</ul>
<p>To avoid - or at least mitigate - these problems, we need to remove these cells at the start of the analysis. This step is commonly referred to as quality control (QC) on the cells. (We will use “library” and “cell” rather interchangeably here, though the distinction will become important when dealing with droplet-based data.) We will demonstrate using a small scRNA-seq dataset from  A. T. L. Lun et al. (<a href="https://osca.bioconductor.org/quality-control.html#ref-lun2017assessing">2017</a>), which is provided with no prior QC so that we can apply our own procedures.</p>

