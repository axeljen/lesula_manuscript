## Custom Python scripts used in this manuscript

### Requirements
	pysam

### allele_balance_filter.py
Masks heterozygous calls with minor allele read support below a specified threshold.

	## Usage example
	python3 allele_balance_filter.py -i input.vcf -o output.vcf  --min-ab 0.25

### variant_stats.py
Script to output some general variant stats, including heterozygosity and homozygosity counts.

	## Usage example
	python3 variant_stats.py input.vcf output.stats

### make_homref_vcf.py
This scripts takes a single-sample vcf file and outputs a new vcf file with all calls set to homozygous reference. Used to estimate the total length of accessible genome in the ROH analysis.

	## Usage example
	python3 make_homref_vcf.py -i input.vcf -o output.homref.vcf


### Conversion from vcf to phylip/fasta formant and translation to amino acids

This was done with scripts from my general bioinfo repo:
https://github.com/axeljen/phylogenomics

After cloning this repo, transcripts can be extracted to phylip format using the vcfToMSA.py script, provided that all coding coordinates of the transcript is listed in a text file with four tab-separated columns:
chrom\tstart\tend\tstrand. 
	## usage example
	python3 vcfToMSA.py -v input_vcf.vcf \
		-R transcript_cds.txt \
		--concat \
		--reference reference_genome.fasta

### convergently fixed amino acid differences in gene transcripts were identified using the fixed_differences.py script
This script also make use of the general bioinfo repo, after cloning this the path to that repo should be accordingly modified in the header of this script.

Then, this script takes an alignment and a sample-pop assignment file for the two groups to check. add -t to first translate the alignment to amino acid sequences. The script will also check for internal stop codons.

	## usage example
	python3 fixed_differences.py -a alignment.phy \
		-p popfile.assignments.txt \
		-t \
		-o output.txt







