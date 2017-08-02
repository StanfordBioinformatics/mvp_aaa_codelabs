<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Stanford University All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

*This codelab was made in collaboration with [Google Genomics](https://github.com/googlegenomics)*

# Part 1: Data Preprocessing and BigQuery Setup

This document details the methods used to prepare genomic data for analysis on Google Cloud.  Here we are working with variant data stored in gVCFs.  gVCFs, unlike standard VCFs, contain reference calls for the entire genome.  Our gVCFs have consequtive reference calls with similar quality grouped into reference blocks.  If the gVCF has individual reference calls on each line, you may want to group together into reference blocks.  You can do this by running [this script](https://github.com/StanfordBioinformatics/googva/blob/master/gvcf-mapper-cl.py) on Google Cloud.  VCFs also work fine with the Google Genomics tools.

## Upload data 
First, we uploaded all the vcfs from our local compute cluster to Google Cloud Storage.  A tool called [gsutil](https://cloud.google.com/storage/docs/gsutil) allows up to interact with Google Cloud.  We'll use this to upload our data.

You will need to first create an account with Google Cloud to use Google Storage.  The path `gs://path/to/google/storage` represents a directory within a bucket in Google Storage.

```r
gsutil cp /path/to/local/data gs://path/to/google/storage
```


## Import into genomics variant set
Once our data is in Google Cloud, we import all the vcfs into a genomics [variant set](https://cloud.google.com/genomics/v1/managing-variants).  This makes out data accessible through the Google Genomics API.

Google Genomics provides a toolset to interact with genomic data with [gcloud](https://cloud.google.com/sdk/gcloud/) (gcloud is installed automatically when you install gsutil).

To get started with the genomics tools simply run the following command.  You will be prompted to install the alpha commands.

```r
gcloud alpha genomics
```

Now we can use the genomics tools provided within gcloud.

First we need to create a new dataset.
```r
gcloud alpha genomics datasets create --name my_dataset
```

gcloud will output something like the following.  We'll need to save the id for later use.
```
Created dataset my_dataset, id: 12406857362375913404
```

Next we'll create a variant set within the dataset.  This is where we'll put our data in the next step.

```r
gcloud alpha genomics variantsets create --dataset-id 12406857362375913404
```
This outputs the following.  Now we'll need to save the variant set id.
```
Created variant set id: 14165412073904006532 "", belonging to dataset id: 12406857362375913404
Created [14165412073904006532].
```

Now we have created a space into which we can copy our variant data.  Copying our VCFs into a variant set formats it in an object oriented manner that can be accessed rapidly by the Genomics API.  


Now we can import our gVCFs into the variant set.  Again we use gcloud.

```r
gcloud alpha genomics variants import --variantset-id 14165412073904006532 --source-uris gs://path/to/google/storage*.vcf 
```

This function submits a job to Google Cloud.  The job id is output.
```
done: false
name: operations/CJ3k0-yGHBCmrYC8BRi1ou-Z5OTE4PEB
```

If the VCFs are very large the import job can take a while.  You can use the following command to monitor the job.

```r
gcloud alpha genomics operations describe operations/CJ3k0-yGHBCmrYC8BRi1ou-Z5OTE4PEB
```

When that command returns `done:true` the import has completed.

## Export to BigQuery

### Genome Call Table
Now that we have our data stored in a variant set we can create a BigQuery table where we can query our data.

```r
gcloud alpha genomics variantsets export 14165412073904006532 my_genome_calls_table --bigquery-dataset my_bq_dataset
```

### Multi Sample Variants Table
For some of the queries we want to run it is much easier to have the call for each sample listed explicitly rather than grouped into reference matching blocks.  To accomplish this we create a table resembling a multi sample VCF with the call for every sample listed at every position where at least one sample had a variant.  To create this table we run a [Cloud Dataflow job](https://github.com/StanfordBioinformatics/codelabs/tree/master/Java/PlatinumGenomes-variant-transformation) to transform our genomic call table into a multi sample variant table. For cohorts of size ~2,000 or larger we recommend
[a more recent version of this pipeline](https://github.com/googlegenomics/codelabs/tree/master/Java/PlatinumGenomes-variant-transformation)
that can be configured to emit a schema optimized for large cohorts.

*NOTE:* The code needs to be compiled before running.  Follow the setup instructions for the variant transformation.

```r
java -cp target/non-variant-segment-transformer-*.jar   \
  com.google.cloud.genomics.examples.TransformNonVariantSegmentData   \
  --project=my-project-name   \
  --stagingLocation=gs://my_bucket/dataflow/staging   \
  --datasetId=14165412073904006532   \
  --allReferences   \
  --hasNonVariantSegments    \
  --outputTable=gbsc-gcp-project-mvp:my_bq_dataset.my_multisample_variants_table \
  --numWorkers=200 \
  --runner=DataflowPipelineRunner
```
--------------------------------------------------------
_Next_: [Part 2: Sample-Level QC](./Sample-Level-QC.md)
