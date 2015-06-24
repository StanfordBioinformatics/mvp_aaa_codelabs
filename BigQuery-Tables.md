# BigQuery Table Overview

This document provides an overview of all the BigQuery tables within the MVP Google Cloud project.

## Table types
Two types of tables are used to store variants: genome_call and multi_sample_variant tables.  Here we provide a descripiton of each of the tables that will be repeated throughout this document.

### genome_call
The genome_call tables contain genomic data in a format similar to [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf). Each row represents a single variant, block of reference calls, or block of low quality calls belonging to a single sample.  

### multi_sample_variants 
The multi_sample_variants tables differ from the genome_calls tables in that each variant is only listed once and the call for every sample is listed in that row, regardless of whether it is a variant, reference, or low quality call.  The motivation for creating this table is to simplify some some queries.  For more information on how this table was created and why please see the [documentation](https://github.com/googlegenomics/codelabs/tree/master/Java/PlatinumGenomes-variant-transformation).

## No Quality Control
These tables contain the raw variant data from HaplotypeCaller with no additional filtering performed.

#### genome_calls_no_qc

#### multi_sample_variants_no_qc

## Sequence Level Quality Control
These tables contain variant data with sequence level filtering.  Variants and reference calls not passing filter criteria are flagged with a -1/-1 genotype. These are intermediate tables used for performing sample and variant level quality control.  

Variants with any flag other than `PASS` in the FILTER column are flagged.

Reference calls not meeting the following criteria are flagged:

  * QUAL >= 30
  * MQ >= 30
  * MQ0 < 4

#### genome_calls_seq_qc

#### multi_sample_variants_seq_qc

## Full Quality Control
The data in these tables have had sample and variant level quality control techniques applied.  Any variants or references not passing any quality control test will have the failed test listed in the QC column.  There are two QC columns, one at the variant info level, and a call.QC column containing quality control metrics specific to a callset.

#### genome_calls_full_qc

#### multi_sample_variants_full_qc

## Other tables

#### patient_info
This table contains patient info for the AAA cohort.
