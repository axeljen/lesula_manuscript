import argparse
from pysam import VariantFile
import sys
parser = argparse.ArgumentParser()

parser.add_argument("--vcf", type=str)
parser.add_argument("--output", type=str)

args = parser.parse_args()

vcf = VariantFile(args.vcf)

# should only be one sample, double check and get sample
if len(list(vcf.header.samples)) > 1:
	print("More than one sample, exiting..")
	sys.exit()
sample = list(vcf.header.samples)[0]

# prep output vcf for homozygous genotypes
vcf_homref = VariantFile(args.output, "w", header = vcf.header)

# loop through all the recs and output filtered and homrefs
for rec in vcf.fetch():
	if not None in rec.samples[sample]['GT']:
		homrec = rec
		homrec.samples[sample]['GT'] = (0,0)
		vcf_homref.write(homrec)
