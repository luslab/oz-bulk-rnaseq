rule star_map:
	input:
		index="/home/camp/ziffo/working/oliver/genomes/index/GRCh38.p12_STAR_index",
		sample="/home/camp/ziffo/working/oliver/projects/airals/fastq_files/{sample}.fastq"
	output:
		"/home/camp/ziffo/working/oliver/projects/airals/alignment_STAR/{sample}.sam"
	shell:
		"STAR --genomeDir {input.index} --readFilesIn {input.sample} --outFileNamePrefix {output}_ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 1"
