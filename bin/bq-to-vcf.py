from __future__ import print_function
import sys
import os
import re
import argparse

CHROM=0
POS=1
ID=2
REF=3
ALT=4
QUAL=5
FILTER=6
INFO=7
CALLS=8

def main(file, sample_list=None, one_based=False):
    if not os.path.exists(file):
        print("Couldn't find input file: %s" % file)
        exit(1)
    samples = {}
    if sample_list is not None:
        samples = samples_from_file(sample_list)
        print_header(samples)
    with open (file, "r") as f:
        for line in f:
            if line.startswith("CHROM"):
                continue
            line = line.rstrip()
            list = line.split(",")
            calls = list[CALLS].split("|")
            if not samples:
                samples = make_sample_dict(calls)
                print_header(samples)
            genotypes = get_genotypes(calls, samples)
            list[ID] = clean_name(list[ID])
            if one_based is True:
                list[POS] = convert_to_one_based(list[POS])
            print_line(list, genotypes, samples)

def clean_name(name):
    ids = name.split("|")
    unique_ids = {}
    for i in ids:
        unique_ids[i] = 1
    clean_ids = ",".join(unique_ids.keys())
    if not name:
        clean_ids = '.'
    return clean_ids

def make_sample_dict(calls):
    count = 0
    samples = {}
    for c in calls:
        if not c:
            continue
        sample = c.split(":")[0]
        samples[count] = sample
        count += 1
    return samples

def samples_from_file(file):
    samples ={}
    count = 0
    with open(file) as f:
        for s in f:
            s = s.rstrip()
            samples[count] = s
            count += 1
    return samples

def get_genotypes(calls, samples):
    genotypes = {}
    for c in calls:
        if not c:
            continue
        sample, genotype = c.split(":")
        genotypes[sample] = genotype
    genotypes = get_no_calls(genotypes, samples)
    genotypes = clean_low_quality(genotypes)
    return genotypes

def get_no_calls(genotypes, samples):
    sample_list = samples.values()
    for s in sample_list:
        if not s in genotypes:
            genotypes[s] = './.'
    return genotypes

def clean_low_quality(genotypes):
    for g in genotypes:
        if genotypes[g] == '-1/-1':
            genotypes[g] = './.'
    return genotypes

def convert_to_one_based(position):
    one_based = int(position) + 1
    return str(one_based)

def print_line(list, genotypes, samples):
    statics = list[0:8]
    genotypes_sorted = []
    count = 0
    while count < len(samples):
        genotypes_sorted.append(genotypes[samples[count]])
        count += 1
    print("\t".join(statics) + "\tGT\t" + "\t".join(genotypes_sorted))

def print_header(samples):
    # This is just the standard header from the VCF specification file.
    stock_header = ['##fileformat=VCFv4.2',
                    '##fileDate=20090805',
                    '##source=myImputationProgramV3.1',
                    '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta',
                    '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
                    '##phasing=partial',
                    '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                    '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
                    '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
                    '##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
                    '##FILTER=<ID=q10,Description="Quality below 10">',
                    '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                    '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">']
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    sample_list = []
    count = 0
    while count < len(samples):
        sample_list.append(samples[count])
        count += 1
    print("\n".join(stock_header))
    print("#" + "\t".join(columns) + "\t" + "\t".join(sample_list))

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Convert csv output from BigQuery to vcf.')

    parser.add_argument("--input", default=None,
                                help="csv file to be converted to vcf")
    parser.add_argument("--sample_list", default=None,
                                help="List of samples to include in the vcf.  If not included samples will be "
                                     "collected from the first variant in the input.  If some samples have missing"
                                     "calls you should provide a sample list.")
    parser.add_argument("--convert_to_one_based", default=False, action='store_true',
                                help="Convert zero based coordinates to one based.")

    options = parser.parse_args()
    if options.input is None:
        print("Exiting, specify input.")
        exit(0)
    return options

if __name__ == "__main__":
    options = parse_command_line()
    main(options.input, options.sample_list, options.convert_to_one_based)