---


---

<blockquote>
<h1 id="align-to-reference-transcriptome">Align to reference transcriptome</h1>
</blockquote>
<p>Download <code>Cell Ranger</code> into bin/  from:<br>
<a href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation">https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation</a></p>
<p>Utilise the CAMP shared 10X reference genome located at: <code>/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0</code></p>
<p>Fastqs are located in: <code>/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/</code></p>
<h1 id="cell-ranger">Cell Ranger</h1>
<p><a href="https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/">https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/</a></p>
<h2 id="mkfastq-converts-bcl-files-to-fastq-files.">mkfastq: converts BCL files to fastq files.</h2>
<h2 id="count-quantifies-single-cell-gene-expression-from-fastq-files">count: quantifies single cell gene expression from fastq files</h2>
<pre class=" language-bash"><code class="prism  language-bash"><span class="token comment"># cp fastq files into cell-line specific directory:</span>
<span class="token function">mkdir</span> /camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled
<span class="token function">mkdir</span> AC2_C2MN
<span class="token function">mkdir</span> ACA_GliaMN
<span class="token function">mkdir</span> AC3_C3MN
<span class="token function">mkdir</span> ACD_CB1D

<span class="token function">cd</span> /camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/

<span class="token comment"># go into each of the fastq subfolders and run:</span>
<span class="token function">ls</span> -1 <span class="token keyword">.</span> <span class="token operator">|</span> <span class="token function">egrep</span> <span class="token string">'TAH421A1_*'</span> <span class="token operator">|</span> <span class="token function">xargs</span> <span class="token function">cp</span> -t /camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/AC2_C2MN
<span class="token function">ls</span> -1 <span class="token keyword">.</span> <span class="token operator">|</span> <span class="token function">egrep</span> <span class="token string">'TAH421A2_*'</span> <span class="token operator">|</span> <span class="token function">xargs</span> <span class="token function">cp</span> -t /camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/ACA_GliaMN
<span class="token function">ls</span> -1 <span class="token keyword">.</span> <span class="token operator">|</span> <span class="token function">egrep</span> <span class="token string">'TAH421A3_*'</span> <span class="token operator">|</span> <span class="token function">xargs</span> <span class="token function">cp</span> -t /camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/AC3_C3MN
<span class="token function">ls</span> -1 <span class="token keyword">.</span> <span class="token operator">|</span> <span class="token function">egrep</span> <span class="token string">'TAH421A4_*'</span> <span class="token operator">|</span> <span class="token function">xargs</span> <span class="token function">cp</span> -t /camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/ACD_CB1D
</code></pre>
<p>Run cellranger count in each of these 4 cell directories separately:</p>
<pre class=" language-bash"><code class="prism  language-bash">sbatch --time<span class="token operator">=</span>48:00:00 --wrap <span class="token string">"cellranger count --id=AC2_C2MN_counts --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/AC2_C2MN --sample=TAH421A1"</span> --job-name<span class="token operator">=</span><span class="token string">"AC2_C2MN_counts"</span> -c 16 --mem-per-cpu<span class="token operator">=</span>7000

sbatch --time<span class="token operator">=</span>48:00:00 --wrap <span class="token string">"cellranger count --id=ACA_GliaMN_counts --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/ACA_GliaMN --sample=TAH421A2"</span> --job-name<span class="token operator">=</span><span class="token string">"ACA_GliaMN_counts"</span> -c 16 --mem-per-cpu<span class="token operator">=</span>7000

sbatch --time<span class="token operator">=</span>48:00:00 --wrap <span class="token string">"cellranger count --id=AC3_C3MN_counts --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/AC3_C3MN --sample=TAH421A3"</span> --job-name<span class="token operator">=</span><span class="token string">"AC3_C3MN_counts"</span> -c 16 --mem-per-cpu<span class="token operator">=</span>7000

sbatch --time<span class="token operator">=</span>48:00:00 --wrap <span class="token string">"cellranger count --id=ACD_CB1D_counts --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/outputs/gandhi-patani/doaa.taha/fastq_pooled/ACD_CB1D --sample=TAH421A4"</span> --job-name<span class="token operator">=</span><span class="token string">"ACD_CB1D_counts"</span> -c 16 --mem-per-cpu<span class="token operator">=</span>7000

/camp/lab/luscomben/inputs/babs-gandhi-patani/doaa.taha/asf/SC19137


/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq
/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq/TAH421A1
/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq/TAH421A2
/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq/TAH421A3
/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190806_K00102_0374_AH5KN3BBXY/fastq/TAH421A4

/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190820_K00102_0382_BH7MGCBBXY/fastq

/camp/stp/babs/outputs/gandhi-patani/doaa.taha/asf/SC19137/190927_K00102_0401_BHF357BBXY/fastq

</code></pre>
<p>–id: this is the name of the new directory<br>
–sample : fastq prefixname.<br>
–transcriptome: reference transcriptome</p>
<p>The output includes the .bam and the .bai index.</p>
<p>View the html summary of reads &amp; mapping QC</p>
<h3 id="droputils">DropUtils</h3>
<p>Alternative to the html view from cellranger count. Helps determine if the droplets contain a cell of just ambient RNA (soup).<br>
Run in R wiht library(DropletUtils)</p>
<pre><code>## Resource to create custom 10x reference
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna

## Determine cell type using AUCc
# https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html
</code></pre>
<h1 id="seurat">Seurat</h1>
<p><a href="https://davetang.org/muse/2017/08/01/getting-started-seurat/">https://davetang.org/muse/2017/08/01/getting-started-seurat/</a></p>

