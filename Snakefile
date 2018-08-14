rule refData_hg19Lengths:
	output:
		"reference_data/hg19.genome"
	shell:
		"curl -s https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > {output}"
