from os.path import join, basename
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

# shell.prefix('PATH=/net/snowwhite/home/beckandy/miniconda3/bin:$PATH')

CHROMOSOMES = list(range(1,23))
ALL_SOMES = CHROMOSOMES[:]
ALL_SOMES.extend(["X","Y"])

ANCESTRAL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2"
ANCESTRALBASE = basename(ANCESTRAL)

REFERENCEDIR = "reference_data"

rule all:
	input:
		"reference_data/gc10kb.bed"

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
		temp("reference_data/human_g1k_v37_mask/human_g1k_v37.premask.fasta")
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
		"tar -xjf {input} -C reference_data"


rule refData_fixedWidthWindows:
	input:
		"reference_data/hg19.genome"
	output:
		win1000 = "reference_data/genome.1000kb.sorted.bed",
		win5000 = "reference_data/genome.5000kb.sorted.bed",
		win100 = "reference_data/genome.100kb.sorted.bed",
		win10 = "reference_data/genome.10kb.sorted.bed",
		winFull = "reference_data/genome.full.sorted.bed"
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
		fasta = "reference_data/human_g1k_v37/human_g1k_v37.fasta",
		bed = "reference_data/genome.10kb.sorted.bed"
	output:
		"reference_data/gc10kb.bed"
	shell:
		"sed s/chr// {input.bed} | bedtools nuc -fi {input.fasta} -bed - > {output}"

rule refData_cpgIslands:
	output:
		"reference_data/cpg_islands_sorted.bed"
	shell:
		"curl -s  http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt | awk \'NR>1\' | sort -k1,1 -k2,2n > {output}"

rule refData_lamin:
	output:
		"reference_data/lamin_B1_LADS2.bed"
	shell:
		"curl -s  \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz\" | gunzip | awk \'NR>1 {{print $2\"\\t\"$3\"\\t\"$4}}\' | bedtools sort -i - > {output}"

rule refData_DNase:
	output:
		"reference_data/DHS.bed"
	shell:
		"curl -s \"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz\" | gunzip | cut -f1-3 | bedtools sort -i - > {output}"

rule refData_RT:
	output:
		"reference_data/lymph_rep_time.txt"
	shell:
		"curl -s \"http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip\" | gunzip > {output}"

rule refData_recomb:
	output:
		"reference_data/recomb_rate.bed"
	shell:
		"""
		curl -s \"http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw\" > \"reference_data/SexAveraged.bw\"
		bigWigToWig \"reference_data/SexAveraged.bw\" \"reference_data/SexAveraged.wig\"
		echo \"CHR\\tSTART\\tEND\\tRATE\" > {output}
		awk \'NR>1\' \"reference_data/SexAveraged.wig\" | cat >> {output}
		"""

rule refData_getHistone:
	output:
		"reference_data/sort.E062-H3K27ac.bed",
		"reference_data/sort.E062-H3K27me3.bed",
		"reference_data/sort.E062-H3K36me3.bed",
		"reference_data/sort.E062-H3K4me1.bed",
		"reference_data/sort.E062-H3K4me3.bed",
		"reference_data/sort.E062-H3K9ac.bed",
		"reference_data/sort.E062-H3K9me3.bed"
	shell:
		"""
		wget -r -nd -P . --accept-regex \'E062\' https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

		gunzip *.broadPeak.gz
		for f in *.broadPeak; do
			mv -- \"$f\" \"reference_data/${{f%.broadPeak}}.bed\"
		done

		for i in reference_data/E062*.bed; do
			outFile=$(echo $i | sed -E  's/(.*)\/(.*)/\\2/')
			bedtools sort -i $i > \"reference_data/sort.$outFile\"
		done
		"""

rule refData_deNovo_goNL:
	input:
		HTTP.remote("https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt", keep_local=True)
	output:
		"DNMs/GoNL_DNMs.txt"
	run:
		shell("mv {input} {output}")

rule refData_deNovo_goldmann:
	input:
		HTTP.remote("https://media.nature.com/original/nature-assets/ng/journal/v48/n8/extref/ng.3597-S3.xlsx", keep_local=True)
	output:
		"DNMs/goldmann_2016_dnms.xlsx"
	run:
		shell("mv {input} {output}")

rule refData_refSeqExons:
	input:
		HTTP.remote("http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed", keep_local=True)
	output:
		"reference_data/GRCh37_RefSeq_sorted.bed"
	run:
		shell("mv {input} {output}")

rule refData_cytobands:
	input:
		HTTP.remote("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz", keep_local=True)
	output:
		cyto="reference_data/cytoBand.txt"
	run:
		shell("gunzip -c {input} > {output}")

rule refData_maskPctPerBand:
	input:
		cyto="reference_data/cytoBand.txt",
		bed="reference_data/testmask2.bed"
	output:
		"reference_data/testcov.bed"
	shell:
		"bedtools coverage -a {input.bed} -b {input.cyto} > {output}"

rule refData_2xHAR:
	input:
		FTP.remote("ftp://ftp.broadinstitute.org/pub/assemblies/mammals/29mammals/2xHARs.bed", keep_local=True)
	output:
		"reference_data/2xHARs.bed"
	run:
		shell("mv {input} {output}")

rule refData_liftover:
	input:
		HTTP.remote("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver", keep_local=True)
	output:
		"reference_data/liftOver"
	run:
		shell("mv {input} {output}")
		shell("chmod +x {output}")

rule refData_overChain:
	input:
		HTTP.remote("http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz", keep_local=True)
	output:
		"reference_data/hg18ToHg19.over.chain.gz"
	run:
		shell("mv {input} {output}")

rule refData_2xHARUnlifted:
	input:
		liftOver="reference_data/liftOver",
		har="reference_data/2xHARs.bed",
		chain="reference_data/hg18ToHg19.over.chain.gz"
	output:
		har="reference_data/2xHARs.hg19.bed",
		bed="reference_data/unlifted.bed"
	run:
		shell("\"{input.liftOver}\" \"{input.har}\" \"{input.chain}\" \"{output.har}\" \"{output.bed}\"")

rule refData_HARSort:
	input:
		"reference_data/2xHARs.hg19.bed"
	output:
		"reference_data/2xHARs.hg19.sort.bed"
	run:
		shell("bedtools sort -i {input} > {output}")

rule refData_AggarwalaVoight:
	input:
		HTTP.remote("https://media.nature.com/original/nature-assets/ng/journal/v48/n4/extref/ng.3511-S2.xlsx", keep_local=True)
	output:
		"reference_data/AV_rates.xlsx"
	run:
		shell("mv {input} {output}")

rule refData_compIndexRef:
	input:
		"reference_data/{genome}/{genome}.fasta"
	output:
		"reference_data/{genome}/chr{chr}.fasta.gz"
	threads: 1
	shell:
		"""
		samtools faidx {input} {wildcards.chr} | bgzip -c > {output}
		samtools faidx {output}
		"""

rule refData_compAllChroms:
	input:
		expand("reference_data/human_g1k_v37/chr{chr}.fasta.gz", chr=CHROMOSOMES)

rule refDat_cleanMask:
	input: "reference_data/human_g1k_v37_mask/human_g1k_v37.premask.fasta"
	output: "reference_data/human_g1k_v37_mask/human_g1k_v37_mask.fasta"
	shell:
		"""
		perl -ane 'if(/\>/){{$a++;print \">$a dna:chromosome\n\"}}else{{print;}}' {input} > {output}
		"""

rule refData_compAllChromsMask:
	input:
		expand("reference_data/human_g1k_v37_mask/chr{chr}.fasta.gz", chr=CHROMOSOMES)

rule refData_compIndexAnc:
	input:
		"reference_data/human_ancestor_GRCh37_e59/human_ancestor_{chr}.fa"
	output:
		"reference_data/human_ancestor_GRCh37_e59/human_ancestor_{chr}.fa.gz"
	threads: 1
	shell:
		"""
		cat {input} | sed sed \"s,^>.*,>$i,\" | bgzip -c > {output}
		samtools faidx {output}
		"""
