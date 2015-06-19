# Data Preparation on Google Cloud

This document describes the process used to create a BigQuery variant table from gVCFs stored in a Google Cloud bucket.

## Converting expanded gVCFs to gVCFs with reference matching blocks

This process was performed using a hadoop streaming job with [gvcf-mapper.py](https://github.com/StanfordBioinformatics/googva/blob/master/gvcf-mapper.py).  This conversion filters each position based on quality and creates references matching blocks for all calls that passed quality and match the reference call.

## Quality Filters

Variants:

  * Filter = PASS

Refererence calls:
  * MQ > 30
  * MQ0 < 4
  * QUALITY > 30

Information about these metrics can be found [here](http://samtools.github.io/hts-specs/VCFv4.1.pdf).

## Commands

#### Set up Hadoop cluster on Google Cloud

```r
$ diff bdutil_env.sh bdutil_env_original.sh
25c25
< CONFIGBUCKET="gbsc-gcp-project-mvp-va_aaa_jobs"
---
  > CONFIGBUCKET=""
28c28
< PROJECT="gbsc-gcp-project-mvp"
---
  > PROJECT=""
50,51c50,51
< # The number of worker nodes in the cluster
  < NUM_WORKERS=40
---
  > # The number of worker nodes in the cluster.
  > NUM_WORKERS=2

./bdutil deploy

gcutil push hadoop-m /Users/gmcinnes/GitHub/googva/gvcf-mapper.py /home/gmcinnes
gcutil push hadoop-m /Users/gmcinnes/GitHub/googva/CustomMultiOutputFormat.java /home/gmcinnes

```

#### Launch Conversion Script 

Now we are ready to log in to the hadoop cluster and run the job.

```r
gcutil ssh hadoop-m

sudo apt-get install -y openjdk-7-jdk screen less vim

javac -cp $(hadoop classpath) -d . CustomMultiOutputFormat.java

jar cvf custom.jar com/custom/CustomMultiOutputFormat.class com/custom/CustomMultiOutputFormat\$LineRecordWriter.class 

sudo su

cd /home/hadoop/hadoop-install

screen

./bin/hadoop jar contrib/streaming/hadoop-streaming-1.2.1.jar \
  -D mapred.max.map.failures.percent=25 \
  -D mapred.max.reduce.failures.percent=25 \
  -D mapred.map.max.attempts=6 \
  -libjars /home/gmcinnes/custom.jar \
  -outputformat com.custom.CustomMultiOutputFormat \
  -input  gs://gbsc-gcp-project-mvp-va_aaa/data/LP6005038-DNA_{A01,A02,A03,B01,B02}/*/vcfs/* \
  -mapper /home/gmcinnes/gvcf-mapper.py \
  -file /home/gmcinnes/gvcf-mapper.py  \
  -reducer org.apache.hadoop.mapred.lib.IdentityReducer \
  -output gs://gbsc-gcp-project-mvp-va_aaa_hadoop/out/5-genome-test/reduced-gvcf-filtered-no-calls-1-based
```

Job information can be found [here](./jobs/gvcf-mapper_5-genomes.html).


## Import To BigQuery

After the gVCFs have been formatted and filtered we can import them into BigQuery.  To do so we need to create a dataset, import the vcf files into the dataset using the Google Genomics API, and finally export the variant dataset to a BigQuery table.

*The following commands were executed from a local terminal.*


#### Create dataset
```r
java -jar genomics-tools-client-java-v1beta2.jar createdataset --project_number 963911152157 --name va-aaa_5-genome-test-dataset
The new dataset was created: {
  "id" : "10441340382221677722",   # We need to save this id for use in the next step
  "isPublic" : false,
  "name" : "va-aaa_5-genome-test-dataset",
  "projectNumber" : "963911152157"
}
```

#### Import variants into dataset

```r
java -jar genomics-tools-client-java-v1beta2.jar importvariants --client_secrets_filename /Users/gmcinnes/client_secrets.json --vcf_file gs://gbsc-gcp-project-mvp-va_aaa_hadoop/out/5-genome-test/reduced-gvcf-filtered/LP*/part* --variant_set_id 10441340382221677722 --poll
``` 

#### Export variant dataset to BigQuery

```r
java -jar genomics-tools-client-java-v1beta2.jar exportvariants --variant_set_id 10441340382221677722 --project_number 963911152157 --bigquery_dataset va_aaa_pilot_data --bigquery_table 5_genome_test_gvcfs --poll

```


## gVCF to VCF

We also want to create a table with only the variant calls.  This allows us to save time and reduce the cost of the queries where we are only interested in variant calls.

#### Export BigQuery Table

The first step is to export the BigQuery gVCF table to a json stored within a bucket in our project.  To do so you click on the arrow next to the table name in the BigQuery console and select 'Export Table'.  A window will pop up titled 'Export to Google Storage'.  

Export format: JSON (Newline Delimited)
Compression: None
Google Cloud Storage URI: gs://gbsc-gcp-project-mvp-va_aaa/5-genome-test-gvcf-BQ-export/*

Instructions for exporting a BigQuery table can be found [here](https://cloud.google.com/bigquery/bigquery-web-ui#exportdata).

Once the data has been exported launch another cluster to perform the analysis.

#### Create a cluster configured for BigQuery
```r
# Create cluster with 40 workers
diff bdutil_env.sh bdutil_env_original.sh
< # The number of worker nodes in the cluster
< NUM_WORKERS=40
---
> # The number of worker nodes in the cluster.
> NUM_WORKERS=2

./bdutil deploy -e bigquery_env.sh

gcutil push hadoop-m /Users/gmcinnes/GitHub/codelabs/Python/PlatinumGenomes-variant-transformation/gvcf* /home/gmcinnes
gcutil push hadoop-m /Users/gmcinnes/GitHub/codelabs/Python/PlatinumGenomes-variant-transformation/platinum_genomes.variants.schema /home/gmcinnes

gcutil ssh hadoop-m

sed -i "s/BIG_QUERY_SOURCE = True/BIG_QUERY_SOURCE = False/" gvcf-expand-mapper.py
sed -i "s/BIG_QUERY_SINK = True/BIG_QUERY_SINK = False/" gvcf-expand-reducer.py

hadoop jar /home/hadoop/hadoop-install/contrib/streaming/hadoop-streaming-1.2.1.jar -input gs://gbsc-gcp-project-mvp-va_aaa/5-genome-test-gvcf-BQ-export/* -file gvcf_expander.py -mapper gvcf-expand-mapper.py -file gvcf-expand-mapper.py  -reducer gvcf-expand-reducer.py -file gvcf-expand-reducer.py -output gs://gbsc-gcp-project-mvp-va_aaa_hadoop/5-genome-test-vcf-expanded

# Create schema for mvp 
bq --project gbsc-gcp-project-mvp show --format json va_aaa_pilot_data.5_genome_test_gvcfs | python -c "import json,sys ; print \"'%s'\" % (json.dumps(json.loads(sys.stdin.readline())['schema']['fields']).replace(\"'\", \"_\"))" > mvp.gvcf.schema

# **The schema needs to be edited before the next step. Remove the quotes from the beginning and end of the file.**

# Export to BQ
bq load --source_format=NEWLINE_DELIMITED_JSON va_aaa_pilot_data.5_genome_test_vcfs gs://gbsc-gcp-project-mvp-va_aaa_hadoop/5-genome-test-vcf-expanded/part* mvp.gvcf.schema

```

Hadoop job information can be found [here](./jobs/Hadoop-job_201502190151_0001.html).




