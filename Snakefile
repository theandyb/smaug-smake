import os 
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

CHROMOSOMES = list(range(1,23))
ALL_SOMES = CHROMOSOMES
ALL_SOMES.extend(["X","Y"])

rule refData_hg19Lengths:
	output:
		"reference_data/hg19.genome"
	shell:
		"curl -s https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > {output}"

rule refData_1000GStrictMask:
	input:
		"reference_data/hg19.genome"
	output:
		"reference_data/testmask2.bed"
	shell:
		"curl -s  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed | "
		"bedtools complement -i - -g {input} | bedtools sort | awk 'match($1, /chr[0-9]+$/) {{print $0}}' > {output}"

rule refData_refGenome:
	output:
		"reference_data/human_g1k_v37/human_g1k_v37.fasta"
	shell:
		"""
		set +e
		curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz | gunzip -q -c > {output}
		exitcode=$?
		if [ $exitcode -eq 1 ]
		then
			exit 1
		else
			exit 0
		fi
		"""
rule refData_mask_v37:
	input:
		fasta = "reference_data/human_g1k_v37/human_g1k_v37.fasta",
		bed = "reference_data/testmask2.bed"
	output:
		"reference_data/human_g1k_v37_mask/human_g1k_v37.premask.fasta"
	shell:
		"bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}"

rule refData_ancestralGenome:
	output:
		expand("reference_data/human_ancestor_GRCh37_e59/human_ancestor_{chromosome}.fa", chromosome=ALL_SOMES),
		expand("reference_data/human_ancestor_GRCh37_e59/human_ancestor_{chromosome}.bed", chromosome=ALL_SOMES),
		"reference_data/human_ancestor_GRCh37_e59/README",
		"reference_data/human_ancestor_GRCh37_e59/README.txt",
		"reference_data/human_ancestor_GRCh37_e59/summary.txt"
	shell:
		"curl -s http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2 > "
		"reference_data/human_ancestor_GRCh37_e59.tar.bz2 && "
		"tar -vjxf reference_data/human_ancestor_GRCh37_e59.tar.bz2"

rule refData_fixedWidthWindows:
	input:
		"reference_data/hg19.genome"
	output:
		win1000= "reference_data/genome.1000kb.sorted.bed",
		win5000 = "reference_data/genome.5000kb.sorted.bed",
		win100 = "reference_data/genome.100kb.sorted.bed",
		win10 = "reference_data/genome.10kb.sorted.bed",
		winFull = "reference_data/genome.full.sorted.bed"
	shell:
		"bedtools makewindows -g {input} -w 1000000 | grep \"-Ev _|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win1000} && "
		"bedtools makewindows -g {input} -w 5000000 | grep \"-Ev _|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win5000} && "
		"bedtools makewindows -g {input} -w 100000 | grep \"-Ev _|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win100} && "
		"bedtools makewindows -g {input} -w 10000 | grep \"-Ev _|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win10} && "
		"bedtools makewindows -g {input} -w 3000000000 | grep \"-Ev _|X|Y|M\" | sort -k 1,1 -k2,2n > {output.winFull}" 
	

rule refData_gcContent:
	input:
		fasta = "reference_data/human_g1k_v37/human_g1k_v37.fasta",
		bed = "reference_data/genome.10kb.sorted.bed"
	output:
		"reference_data/gc10kb.bed"
	shell:
		"sed s/chr// {input.bed} | bedtools nuc -fi {input.fasta} -bed - > {output}"
