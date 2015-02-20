# Data Preparation on HPC Cluster

This document describes the process used to prepare variant data for analysis on an HPC cluster.  The cluster used in the scenario utilized a Sun Grid Engine job queueing system.  

In this process we wish to prepare data from five genomes for analysis.  The initial VCFs were split by chromosome for each genome.  We want to create a multisample VCF with a column for each genome so we can run some analysis tools that look at statistics accross genomes.  Once a multisample VCF is generated all quality information is discarded except for the quality information corresponding the the first genome.  It is important that we do not include low quality positions in our analysis, so we must first filter each VCF individually.  

## VCF Filtering

Each VCF needed to be filtered to meet quality criteria.  Variants and reference calls were processed with different criteria.  This requires the VCFs to be split into a VCF with only reference calls and another with only variant calls.  This was done with done with a bash script that runs [vcf-filter](https://github.com/ekg/vcflib) on the reference calls and [vcftools](http://vcftools.sourceforge.net/) for the variant calls.  A [bash script](./bin/filter_vcfs.sh) was used to run the filtering and merge the output together, creating a single filtered VCF.  The resulting VCF is compressed and indexed using [bgzip and tabix](http://samtools.sourceforge.net/tabix.shtml).

#### Filtering Cutoffs

Variants:

  * FILTER = PASS

Reference Calls:

  * MQ > 30  
  * MQ0 < 4
  * QUAL > 30

An additional cutoff was added for the reference calls due to an artifact in the VCFs.  In some cases when a variant is called there is a '.' in the GENOTYPE field.  

  * AC = 0

Information about these metrics can be found [here](http://samtools.github.io/hts-specs/VCFv4.1.pdf).

```r
# Submitting the jobs
# For each genome in our sample set and for each chromosome, submit a job to the cluster to filter the corresponding VCF.
for s in LP6005038-DNA_A01 LP6005038-DNA_B02 LP6005038-DNA_A03 LP6005038-DNA_B01 LP6005038-DNA_A02 ; do for chr in $(seq 1 22) Y M; do qsub -cwd -b y -N filter_vcf -o $s.chr$chr.out -e $s.chr$chr.err bash ~/bin/filter_vcfs.sh $s /srv/gsfs0/projects/mvp/Pilot/reproc-out/$s/*/vcfs/chr$chr.vcf.gz ; done ; done

#Example log
cat LP6005038-DNA_A01.chr1.out
Sample: LP6005038-DNA_A01
Chr: chr1
Filtering variants...
Filtering reference calls...
Merging vcfs...
Cleaning up...
Finished.
```

## Creating a multisample VCF

#### Merging by chromosome

In order to run analysis across many genomes we need to create a multi sample VCF.  We use [bcftools](http://samtools.github.io/bcftools/bcftools.html) to merge the VCFs by chromosome.

```r
# Example input file
cat chr1.sample_set.txt
../filtered/LP6005038-DNA_A01/LP6005038-DNA_A01.chr1.filtered.vcf.gz
../filtered/LP6005038-DNA_A02/LP6005038-DNA_A02.chr1.filtered.vcf.gz
../filtered/LP6005038-DNA_A03/LP6005038-DNA_A03.chr1.filtered.vcf.gz
../filtered/LP6005038-DNA_B01/LP6005038-DNA_B01.chr1.filtered.vcf.gz
../filtered/LP6005038-DNA_B02/LP6005038-DNA_B02.chr1.filtered.vcf.gz

# Submitting the jobs
# For each input file (one for each chromosome) submit a job to merge all the VCFs listed within the file into a multi sample VCF.
for file in chr*.sample_set.txt ; do chr=`echo $file | cut -f1 -d.` ; qsub -cwd -b y -l h_vmem=5G /srv/gs1/software/bcftools/1.2/bin/bcftools merge -O z -o $chr.merged.sample_set.vcf.gz -l $file ; done
```

#### Concatenating chromsome VCFs into a genome VCF

Now we have a multi sample VCF for each chromsome with a column for each of our five genomes.  Now we will concatenate each of these chomosome VCFs into a single multisample VCF for the entire genome.  We use vcftools for this task.

```r
qsub -N vcf-concat -o chrALL.merged.sample_set.vcf -l h_vmem=10G -cwd -b y "module load vcftools ; vcf-concat -f chr_merged.txt"
```

#### BRCA1 merged VCF

We also want a merged VCF specifically for BRCA1.  Again, we use vcftools.

```r
vcftools --chr chr17 --from-bp 41196312 --to-bp 41277500 --recode --gzvcf chr17.merged.sample_set.vcf.gz --stdout | bgzip > brca1.merged.sample_set.vcf.gz

tabix brca1.merged.sample_set.vcf.gz
```