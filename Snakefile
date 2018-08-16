from os.path import join, basename
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

CHROMOSOMES = list(range(1,23))
ALL_SOMES = CHROMOSOMES
ALL_SOMES.extend(["X","Y"])

ANCESTRAL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2"
ANCESTRALBASE = basename(ANCESTRAL)


rule all:
	input:
		"reference_dir/gc10kb.bed"

rule refData_hg19Lengths:
	output:
		"reference_dir/hg19.genome"
	shell:
		"curl -s https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > {output}"

rule refData_1000GStrictMask:
	input:
		"reference_dir/hg19.genome"
	output:
		"reference_dir/testmask2.bed"
	shell:
		"curl -s  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed | "
		"bedtools complement -i - -g {input} | bedtools sort | awk 'match($1, /chr[0-9]+$/) {{print $0}}' > {output}"

rule refData_refGenome:
	output:
		"reference_dir/human_g1k_v37/human_g1k_v37.fasta"
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
		fasta = "reference_dir/human_g1k_v37/human_g1k_v37.fasta",
		bed = "reference_dir/testmask2.bed"
	output:
		"reference_dir/human_g1k_v37_mask/human_g1k_v37.premask.fasta"
	shell:
		"bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output}"

rule refData_getAncestral:
	input:
		FTP.remote(ANCESTRAL, keep_local=True)
	output:
		temp(join(REFERENCEDIR, ANCESTRALBASE))
	shell:
		"mv {input} {output}"

rule refData_decompressAncestral:
	input:
		join(REFERENCEDIR, ANCESTRALBASE)
	output:
		directory(join(REFERENCEDIR, ANCESTRALBASE.replace(".tar.bz2", "")))
	shell:
		"tar -xjf {input} -C reference_dir"


rule refData_fixedWidthWindows:
	input:
		"reference_dir/hg19.genome"
	output:
		win1000 = "reference_dir/genome.1000kb.sorted.bed",
		win5000 = "reference_dir/genome.5000kb.sorted.bed",
		win100 = "reference_dir/genome.100kb.sorted.bed",
		win10 = "reference_dir/genome.10kb.sorted.bed",
		winFull = "reference_dir/genome.full.sorted.bed"
	shell:
		"""
		bedtools makewindows -g {input} -w 1000000 | grep -Ev \"_|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win1000}
		bedtools makewindows -g {input} -w 5000000 | grep -Ev \"_|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win5000}
		bedtools makewindows -g {input} -w 100000 | grep -Ev \"_|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win100}
		bedtools makewindows -g {input} -w 10000 | grep -Ev \"_|X|Y|M\" | sort -k 1,1 -k2,2n > {output.win10}
		bedtools makewindows -g {input} -w 3000000000 | grep -Ev \"_|X|Y|M\" | sort -k 1,1 -k2,2n > {output.winFull}
		"""

rule refData_gcContent:
	input:
		fasta = "reference_dir/human_g1k_v37/human_g1k_v37.fasta",
		bed = "reference_dir/genome.10kb.sorted.bed"
	output:
		"reference_dir/gc10kb.bed"
	shell:
		"sed s/chr// {input.bed} | bedtools nuc -fi {input.fasta} -bed - > {output}"

rule refData_cpgIslands:
	output:
		"reference_dir/cpg_islands_sorted.bed"
	shell:
		"curl -s  http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt | awk \'NR>1\' | sort -k1,1 -k2,2n > {output}"

rule refData_lamin:
	output:
		"reference_dir/lamin_B1_LADS2.bed"
	shell:
		"curl -s  \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz\" | gunzip | awk \'NR>1 {print $2\"\t\"$3\"\t\"$4}\' | bedtools sort -i - > lamin_B1_LADS2.bed"
