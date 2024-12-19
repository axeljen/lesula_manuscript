import argparse
from pysam import VariantFile

parser = argparse.ArgumentParser(description="Mask heterozygous calls with skewed allele balance in a vcf file.")

parser.add_argument('-i', '--input', type=str, help='Input vcf', required=True)
parser.add_argument('-o', '--output', type=str, help='Output vcf (can be gzipped).', required=True)

#set heterozygotes with minor allele support below specified threashold
parser.add_argument('--min-ab', type=float, help="Minimum support for minor allele in heterozygotes.")

args = parser.parse_args()


def min_AB_filter(input_vcf, output_vcf, AB_threashold):
    with VariantFile(input_vcf) as vcf:
        output = VariantFile(output_vcf, 'w', header=vcf.header)
        samples = list(vcf.header.samples)
        for rec in vcf.fetch():
            if len(rec.alleles) > 1:
                for sample in samples:
                    if rec.samples[sample]['GT'][0] != rec.samples[sample]['GT'][1]:
                        # if heterozygous, get reads supporting both alleles
                        AB_ref,AB_alt = rec.samples[sample]['AD'][0], rec.samples[sample]['AD'][1]
                        if AB_ref <= AB_alt:
                            try:
                                # if ref allele is the minor, get allele support for this as minor support
                                AB_minor = AB_ref / (AB_ref + AB_alt)
                            except:
                                # occasionally fails if there are no reads for the ref allele, set these to missing
                                rec.samples[sample]['GT'] = (None,None)
                        else:
                            # otherwise get it for the alternate allele
                            AB_minor = AB_alt / (AB_ref + AB_alt)
                        if not AB_minor >= AB_threashold:
                            # if minor allele support doesn't pass, set genotype to none
                            rec.samples[sample]['GT'] = (None,None)
                output.write(rec)

# run the ab-filtration
if args.min_ab:
    min_AB_filter(args.input, args.output, args.min_ab)