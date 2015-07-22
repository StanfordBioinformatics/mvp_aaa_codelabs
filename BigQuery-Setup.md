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

This document details the methods used to prepare genomic data for analysis on Google Cloud

## Upload data 
First, we uploaded all the vcfs from our local compute cluster to Google Cloud Storage.  A [bash script](./bin/gc_upoad.sh) was used to to execute the upload for each sample.  

```r
./gc_upload.sh sample_list
```

See [here](https://cloud.google.com/storage/docs/gsutil) for more information on gsutil.

## Sequence level filtering and gvcf conversion
Our variant data for this project is in gVCF format with every genomic position listed.  We want to compress all sequenctial reference calls into reference matching blocks in order to save time and compute cost when querying the data in BigQuery later.  To do this we use a [MapReduce script](https://github.com/StanfordBioinformatics/googva/blob/master/gvcf-mapper.py) running on Hadoop streaming.  This script also flags low quality reference calls and variants.  Flagged reference calls and variants are assigned a genotype of -1/-1 and sequential low quality calls are grouped in to low quality blocks.

Variants are requried to have `PASS` in the FILTER column.  Variants with any other descriptors are flagged.

Referece calls are flagged based on the following metrics:

  * QUAL >= 30
  * MQ >= 30
  * MQ0 < 4
  
First, we need to set up a compute cluster on Google Cloud.  For information on using bdutil see [here](https://cloud.google.com/hadoop/bdutil). For this job we'll use 100 worker nodes, so we set `NUM_WORKERS` in bdutil_env.sh to 100.
```r
./bdutil deploy
```

And copy the files we'll need.
```r
gcloud compute copy-files /Users/gmcinnes/src/googva/gvcf-mapper.py hadoop-m:~/
gcloud compute copy-files /Users/gmcinnes/src/googva/CustomMultiOutputFormat.java hadoop-m:~/
```

Now we can ssh into the cluster and set things up from there.
```r
gcloud compute ssh hadoop-m

sudo apt-get install -y openjdk-7-jdk screen less vim

javac -cp $(hadoop classpath) -d . CustomMultiOutputFormat.java

jar cvf custom.jar com/custom/CustomMultiOutputFormat.class com/custom/CustomMultiOutputFormat\$LineRecordWriter.class 

sudo su

cd /home/hadoop/hadoop-install

screen

./bin/hadoop jar contrib/streaming/hadoop-streaming-1.2.1.jar \
  -libjars /home/gmcinnes/custom.jar \
  -outputformat com.custom.CustomMultiOutputFormat \
  -input  gs://gbsc-gcp-project-mvp-va_aaa/data/LP*/*/vcfs/* \
  -mapper /home/gmcinnes/gvcf-mapper.py \
  -file /home/gmcinnes/gvcf-mapper.py  \
  -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
  -output gs://gbsc-gcp-project-mvp-va_aaa_hadoop/out/reduced-gvcfs-460-genomes
```

## Import into genomics variant set
Next, we import all the vcfs into a genomics [variant set](https://cloud.google.com/genomics/v1beta2/managing-variants).  This makes out data accessible through the genomics API.

First we need to create a new variant set.
```r
java -jar genomics-tools-client-java-v1beta2.jar createdataset \
  --project_number 963911152157 \
  --name AAA_460_Genomes
```
```
The new dataset was created: {
  "id" : "12721125545898647337",
  "name" : "AAA_460_Genomes",
  "projectNumber" : "963911152157"
}
```

We need to save the id of the new dataset (12721125545898647337) for use in future steps.

We then import the vcfs into the new dataset.

```r
java -jar genomics-tools-client-java-v1beta2.jar importvariants \
  --variant_set_id 12721125545898647337 \
  --vcf_file gs://gbsc-gcp-project-mvp-va_aaa_hadoop/out/reduced-gvcfs-460-genomes/LP*/*
```

## Export to BigQuery

### Genome Call Table
Now that we have our data stored in a variant set we can create a BigQuery table where we can query our data.  We name the table genome_calls_no_qc to indicate that this table has not yet undergone quality control.  We will go over this in the next step.

```r
java -jar genomics-tools-client-java-v1beta2.jar exportvariants \
  --variant_set_id 12721125545898647337 \
  --project_number 963911152157 \
  --bigquery_dataset va_aaa_genomic_data \ 
  --bigquery_table genome_calls_seq_qc
```

### Multi Sample Variants Table
For some of the queries we want to run it is much easier to have the call for each sample listed explicitly rather than grouped into reference matching blocks.  To accomplish this we create a table resembling a multi sample VCF with the call for every sample listed at every position where at least one sample had a variant.  To create this table we run a [Cloud Dataflow job](https://github.com/StanfordBioinformatics/codelabs/tree/master/Java/PlatinumGenomes-variant-transformation) to transform our genomic call table into a multi sample variant table.

*NOTE:* The code needs to be compiled before running.  Follow the setup instructions for the variant transformation.

```r
java -cp target/non-variant-segment-transformer-*.jar   \
  com.google.cloud.genomics.examples.TransformNonVariantSegmentData   \
  --project=gbsc-gcp-project-mvp   \
  --stagingLocation=gs://gbsc-gcp-project-mvp-va_aaa_hadoop/dataflow/staging   \
  --genomicsSecretsFile=/Users/gmcinnes/bin/google_genomics/client_secrets.json   \
  --datasetId=12721125545898647337   \
  --allReferences   \
  --hasNonVariantSegments    \
  --outputTable=gbsc-gcp-project-mvp:va_aaa_pilot_data.multi_sample_variants_seq_qc \
  --numWorkers=200 \
  --runner=DataflowPipelineRunner
```
--------------------------------------------------------
_Next_: [Part 2: Sample-Level QC](./Sample-Level-QC.md)
