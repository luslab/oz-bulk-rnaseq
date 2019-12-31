---


---

<blockquote>
<h1 id="overview">Overview</h1>
</blockquote>
<p><img src="https://osca.bioconductor.org/images/Workflow.png" alt="enter image description here"></p>
<p><img src="https://media.springernature.com/full/springer-static/image/art:10.1038/s12276-018-0071-8/MediaObjects/12276_2018_71_Fig1_HTML.jpg?as=webp" alt="enter image description here"></p>
<p><strong>a</strong> The limiting dilution method isolates individual cells, leveraging the statistical distribution of diluted cells.<br>
<strong>b</strong> Micromanipulation involves collecting single cells using microscope-guided capillary pipettes.<br>
<strong>c</strong> FACS isolates highly purified single cells by tagging cells with fluorescent marker proteins.<br>
<strong>d</strong> Laser capture microdissection (LCM) utilizes a laser system aided by a computer system to isolate cells from solid samples.<br>
<strong>e</strong> Microfluidic technology for single-cell isolation requires nanoliter-sized volumes. An example of in-house microdroplet-based microfluidics (e.g., Drop-Seq).<br>
<strong>f</strong> The CellSearch system enumerates CTCs from patient blood samples by using a magnet conjugated with CTC binding antibodies.<br>
<strong>g</strong> A schematic example of droplet-based library generation. Libraries for scRNA-seq are typically generated via cell lysis, reverse transcription into first-strand cDNA using uniquely barcoded beads, second-strand synthesis, and cDNA amplification</p>
<h1 id="library-preparation-methods">Library Preparation Methods</h1>
<p><img src="https://lh3.googleusercontent.com/3qamWTEaXpJquBxrnrJ5lCPiMpRzLg2_IP_bcLOYSYFwA2TS-YYETlQb1H0WvxdvtKmDd7tYQIA7-g" alt="enter image description here"></p>
<h1 id="computational-analysis">Computational Analysis</h1>
<p><a href="https://osca.bioconductor.org/overview.html">https://osca.bioconductor.org/overview.html</a></p>
<p><img src="https://media.springernature.com/full/springer-static/image/art:10.1038/s12276-018-0071-8/MediaObjects/12276_2018_71_Fig2_HTML.jpg?as=webp" alt="enter image description here"></p>
<ol>
<li><strong>Quality control</strong></li>
</ol>
<p>Remove:<br>
low-quality cells which have been damaged during processing or may not have been fully captured by sequencing.<br>
low quality bases<br>
adapter sequences</p>
<p><strong>Metrics:</strong><br>
total counts per cell<br>
proportion of spike-in<br>
mitochondrial reads<br>
number of detected features</p>
<p><strong>Tools:</strong><br>
FastQC <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">http://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></p>
<ol start="2">
<li><strong>Read Alignment</strong></li>
</ol>
<p>Tools:<br>
BWA<br>
STAR</p>
<p>Metrics:<br>
uniquely mapped reads<br>
reads mapped to annotated exons<br>
coverage patterns</p>
<ol start="3">
<li><strong>Nomalisation</strong>:</li>
</ol>
<p>Convert counts into normalized expression values to eliminate cell-specific biases (e.g., in capture efficiency). This allows explicit comparisons across cells. Gene read counts in each cell is proportional to gene-specific expression level &amp; cell-specific scaling factors.</p>
<p>Raw counts are normalised using <strong>scaling factor estimates</strong> by standardising across cells (assuming most genes are not differentially expressed). e.g. RPKM, TPM</p>
<ol start="4">
<li><strong>Filtering aka Feature selection</strong>:</li>
</ol>
<p>Pick a subset of interesting features (genes) for downstream analysis by modelling the variance across cells for each gene and <strong>retaining genes that are highly variable</strong>. The aim is to reduce computational overhead and noise from uninteresting genes.</p>
<ol start="5">
<li><strong>Dimensionality reduction</strong>:  compact the data and further reduce noise. <strong>PCA</strong> is typically used to obtain an initial low-rank representation for more computational work, followed by more aggressive methods like <strong>tSNE</strong> (t-stochastic neighbor embedding) and <strong>UMAP</strong> for visualization.</li>
</ol>
<p><img src="https://media.springernature.com/full/springer-static/image/art:10.1038/s12276-018-0071-8/MediaObjects/12276_2018_71_Fig5_HTML.jpg?as=webp" alt="enter image description here"></p>
<ol start="6">
<li>
<p><strong>Cluster</strong> cells into groups according to similarities in their (normalized) expression profiles. This aims to obtain groupings that serve as empirical proxies for distinct biological states. We interpret these groupings by identifying differentially expressed marker genes between clusters.</p>
</li>
<li>
<p><strong>Cell type identification</strong></p>
</li>
</ol>
<p><img src="https://media.springernature.com/full/springer-static/image/art:10.1038/s12276-018-0071-8/MediaObjects/12276_2018_71_Fig4_HTML.jpg?as=webp" alt="enter image description here"><strong>a</strong> Technical batch effects are a well-known problem in scRNA-seq when the experiment (condition) is conducted in different plates (environment). Cell-specific scaling factors, such as capture and RT efficiency, dropout/amplification bias, dilution factor, and sequencing amount, must be considered in the normalization step.<br>
<strong>b</strong> Single-cell latent variable model (scLVM) can effectively remove the variation explained by the cell-cycle effect. The clear separation is lost in scLVM-corrected expression data using PCA (visualization adapted from ref. <a href="https://www.nature.com/articles/s12276-018-0071-8#ref-CR58" title="Buettner, F. et al. Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nat. Biotechnol. 33, 155–160 (2015).">58</a>).<br>
<strong>c</strong> The expression value <em>y</em> can be modelled as a linear combination of <em>r</em> technical and biological factors and <em>k</em> latent factors with a noise matrix</p>

