import argparse
import sys
## This path should be modified to point to the location of the phylogenomics repo
sys.path.append('/path/to/phylogenomics')
import functions as fn

# function to convert a fasta sequence to an amino acid sequence
def translate(seq, break_on_stop = False):
	# make a list of codons by splitting the sequence in triplets
	codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
	aaseq = ''
	# loop over the codons and translate them
	for codon in codons:
		aa = fn.codon2aa(codon)
		if aa == "*" and break_on_stop:
			break
		aaseq += aa
	return aaseq


# function to check if the two groups are fixed for a change
def checkFixedChange(pops,aln, pos,min_samples):
	# fetch all alleles at this position for both groups
	alleles = {}
	missingness = {p: 0 for p in pops}
	for pop in pops:
		if not pop in alleles.keys():
			alleles[pop] = []
		for sample in pops[pop]:
			allele = aln.sequences[sample].sequence[pos]
			if not allele in ['-', '?', 'X']:
				alleles[pop].append(allele)
			else:
				missingness[pop] += 1
		# if either pop only has missing data, return false
		if len(pops[pop]) - missingness[pop] < min_samples:
			return False,'missing', alleles
	for allele in alleles[list(pops.keys())[0]]:
		# if the allele is in the other group, return false
		if allele in alleles[list(pops.keys())[1]]:
			return False, 'shared', alleles
	# if we get here, the two groups are fixed for a change
	print(alleles)
	return True, 'fixed', alleles

# missing data per sample
def missingSites(aln, sample):
	missing = 0
	for pos in range(aln.length):
		if aln.sequences[sample].sequence[pos] in ['-', '?', 'X','n']:
			missing += 1
	return missing

# prep a parser for the command line arguments
parser = argparse.ArgumentParser(description='Identify/count fixed amino acid changes between two groups, based on either an aa or nt sequence.')

# popfile with samples and groups for the samples to consider
parser.add_argument('-p', '--popfile', type=str, required=True, help='Popfile with samples and groups for the samples to consider.')
# fasta file with sequences
parser.add_argument('-a', '--alignment', type=str, required=True, help='Alignment file with sequences, or list of files.')
# output file
parser.add_argument('-o', '--output', type=str, required=True, help='Output file.')
# store true argument to check if we should translate or not
parser.add_argument('-t', '--translate', action='store_true', help='Translate the sequences to amino acids before comparing them.')
# minimum number of remaining sequences after filtration per group
parser.add_argument('-m', '--min-samples', type=int, help='Minimum number of remaining sequences after filtration per group.', default = 1)
# arguments specifying the groups to contrast against each other
parser.add_argument('-p2', '--pop2', type=str, help="Name of pop1 as labelled in popfile.", default = None)
parser.add_argument('-p1', '--pop1', type=str, help="Name of pop2 as labelled in popfile.", default = None)

parser.add_argument('--min-sites', type=int, help="Minimum number of good sites for output.", default = 10)


# parse the arguments
args = parser.parse_args()

if args.alignment.endswith(("fa","fasta","phy","phylip")):
	alignments = [args.alignment]
else:
	with open(args.alignment) as f:
		alignments = [line.strip() for line in f]

# read the popfile
pops = fn.parsePopfile(args.popfile)
# and make a list of all the samples to keep from this
samples = []
for pop in pops:
	samples += pops[pop]

results = []

for alignment in alignments:
	# read the alignment 
	aln = fn.readSequenceFile(alignment)

	# prune out the samples that are not in the popfile
	aln.subsetSamples(samples)

	# if translate is true, do this now
	if args.translate:
		aln.nt2aa()

	# check if any sequences contain premature stop codons
	premature_stop_codon = {pop: 0 for pop in pops}
	for pop in pops:
		for sample in pops[pop]:
			# check if the sequence contains a stop codon (before the last position)
			if '*' in aln.sequences[sample].sequence[:-1]:
				premature_stop_codon[pop] += 1
				print('Sample {} in group {} contains a premature stop codon.'.format(sample, pop))
	# convert stop codon counts to frequencies
	premature_stop_codon_freq = {pop: premature_stop_codon[pop] / len(pops[pop]) for pop in pops}

	# make a dictionary to store the fixed changes
	fixed_changes = {}
	fixed_changes_count = 0
	missing_data = 0
	evaluated_sites = 0

	# loop through all the positions in the alignment
	for pos in range(aln.length):
		# check if the two groups are fixed for a change
		fixed, reason, alleles = checkFixedChange(pops, aln, pos, min_samples=args.min_samples)
		# if they are, add this to the dictionary
		if fixed:
			fixed_changes[pos] = (reason, alleles)
			fixed_changes_count += 1
			evaluated_sites += 1
		elif reason == 'missing':
			evaluated_sites += 1
			missing_data += 1
		elif reason == 'shared':
			evaluated_sites += 1

	# print the results
	print('Evaluated sites: %i' % evaluated_sites)
	print('Fixed changes: %i' % fixed_changes_count)
	print('Missing data: %i' % missing_data)

	# get length of the alignment
	aln_len = aln.length
	# check each sample for missing data
	missingness_per_sample = {}
	print("Checking sample missingness.")
	for sample in aln.sequences:
		missing = missingSites(aln, sample) / aln_len
		missingness_per_sample[sample] = missing
	print("Done, writing output.")

	# if any population are fixed for internal stop codons, don't write any fixed site output
	for pop in premature_stop_codon:
		if premature_stop_codon[pop] / len(pops[pop]) == 1:
			fixed_changes = {}
	results.append({'alignment': alignment, 'pop1': args.pop1, 'pop2': args.pop2, 'fixed_changes': fixed_changes_count, 'missing_data': missing_data, 'evaluated_sites': evaluated_sites, 'p1_premature_stops': premature_stop_codon[args.pop1], 'p2_premature_stops': premature_stop_codon[args.pop2]})


with open(args.output, 'w') as f:
	f.write("file\tp1\tp2\tfixed_changes\tmissing_data\tevaluated_sites\tp1_premature_stops\tp2_premature_stops\n")
	for result in results:
		f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(result['alignment'], result['pop1'], result['pop2'], result['fixed_changes'], result['missing_data'], result['evaluated_sites'], result['p1_premature_stops'], result['p2_premature_stops']))
