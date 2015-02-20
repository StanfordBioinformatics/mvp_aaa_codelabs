#!/bin/bash

set -e 

# Load modules
module load vcftools;
module load tabix;
module load vcflib;

# Set input
sample=$1;
file=$2;
chr=`echo ${file##*/} | cut -f1 -d.`
output_path=/srv/gsfs0/projects/mvp/Pilot/merged_vcfs/filtered/$sample
var_vcf=$output_path/$sample.$chr.variants.vcf.gz
ref_vcf=$output_path/$sample.$chr.refs.vcf.gz
filtered_vcf=$output_path/$sample.$chr.filtered.vcf.gz

# Check if directory exists
if [ ! -d $output_path ]; then
	mkdir $output_path
fi

echo "Sample: $sample"
echo "Chr: $chr"

# Filter variants
echo "Filtering variants..."
vcftools --gzvcf $file --keep-filtered PASS --recode --stdout | bgzip > $var_vcf
tabix $var_vcf

# Filter reference matching sites
echo "Filtering reference calls..."
vcffilter -f "AC = 0 & MQ > 30 & QUAL > 30 & MQ0 < 4" -g "GT = 0/0"  $file | bgzip > $ref_vcf
tabix $ref_vcf

# Merge into single vcf
echo "Merging vcfs..."
vcf-concat $var_vcf $ref_vcf | vcf-sort | bgzip > $filtered_vcf
tabix $filtered_vcf

# Remove intermediate files
echo "Cleaning up..."
rm $var_vcf
rm $var_vcf.tbi
rm $ref_vcf
rm $ref_vcf.tbi

echo "Finished."
