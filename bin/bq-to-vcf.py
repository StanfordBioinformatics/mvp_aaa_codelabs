from __future__ import print_function
import sys
import os
import re


CHROM=0
POS=1
ID=2
REF=3
ALT=4
QUAL=5
FILTER=6
INFO=7
CALLS=8

def main():
    file = '/Users/gmcinnes/data/sample-bq-vcf.csv'
    samples = {}
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
            genotypes = get_genotypes(calls)
            list[ID] = clean_name(list[ID])
            print_line(list, genotypes, samples)

def clean_name(name):
    ids = name.split("|")
    unique_ids = {}
    for i in ids:
        unique_ids[i] = 1
    clean_ids = ",".join(unique_ids.keys())
    return clean_ids

def make_sample_dict(calls):
    count = 0
    samples = {}
    for c in calls:
        sample = c.split(":")[0]
        samples[count] = sample
        count += 1
    return samples

def get_genotypes(calls):
    genotypes = {}
    for c in calls:
        sample, genotype = c.split(":")
        genotypes[sample] = genotype
    return genotypes

def print_line(list, genotypes, samples):
    statics = list[0:8]
    genotypes_sorted = []
    count = 0
    while count < len(samples):
        genotypes_sorted.append(genotypes[samples[count]])
        count += 1
    print("\t".join(statics) + "\tGT\t" + "\t".join(genotypes_sorted))

def print_header(samples):
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    sample_list = []
    count = 0
    while count < len(samples):
        sample_list.append(samples[count])
        count += 1
    print("#" + "\t".join(columns) + "\t" + "\t".join(sample_list))

if __name__ == "__main__":
    main()