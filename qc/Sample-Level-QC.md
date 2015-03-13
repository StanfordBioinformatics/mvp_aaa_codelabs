<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

# Part 3: Sample-Level QC





In Part 3 of the codelab, we perform some quality control analyses that could help to identify any problematic genomes that should be removed from the cohort before proceeding with further analysis.  The appropriate cut off thresholds will depend upon the input dataset and/or other factors.

* [Missingness Rate](#missingness-rate)
* [Singleton Rate](#singleton-rate)
* [Heterozygosity Rate and Inbreeding Coefficient](#homozygosity-rate-and-inbreeding-coefficient)
* [Sex Inference](#sex-inference)
* [Ethnicity Inference](#ethnicity-inference)
* [Genome Similarity](#genome-similarity)

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)
```

## Missingness Rate

For each genome, determine the percentage of sites explicitly called as a no-call.  If this percentage is too high, the genome may be problematic.


```r
result <- DisplayAndDispatchQuery("./sql/sample-level-missingness.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Compute the ratio of positions corresponding to no-calls versus all positions
# called (reference, variant, and no-calls).
SELECT
  call.call_set_name,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM (
  SELECT
    call.call_set_name,
    SUM(IF(has_no_calls, delta, 0)) AS no_calls,
    SUM(delta) AS all_calls
  FROM (
    SELECT
      END - start AS delta,
      call.call_set_name,
      SOME(call.genotype == -1) WITHIN call AS has_no_calls,
    FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    )
  GROUP BY
    call.call_set_name)
ORDER BY
  call.call_set_name
```
Number of rows returned by this query: 5.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:43:08 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 81237179 </td> <td align="right"> 2147483647 </td> <td align="right"> 0.03 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 97410800 </td> <td align="right"> 2147483647 </td> <td align="right"> 0.03 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 88412773 </td> <td align="right"> 2147483647 </td> <td align="right"> 0.03 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 88709511 </td> <td align="right"> 2147483647 </td> <td align="right"> 0.03 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 83989140 </td> <td align="right"> 2147483647 </td> <td align="right"> 0.03 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=missingness_rate)) +
  theme(axis.text.x=if(nrow(result) <= 20)
    {element_text(angle = 90, hjust = 1)} else {element_blank()}) +
  xlab("Sample") +
  ylab("Missingness Rate") +
  ggtitle("Genome-Specific Missingness")
```

<img src="figure/sampleMissingness-1.png" title="plot of chunk sampleMissingness" alt="plot of chunk sampleMissingness" style="display: block; margin: auto;" />

## Singleton Rate

For each genome, count the number of variants shared by no other member of the cohort.  Too many private calls for a particular individual may indicate a problem.


```r
result <- DisplayAndDispatchQuery("./sql/private-variants.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Compute private variants counts for each sample.
SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS private_variant_count,
FROM (
  SELECT
    reference_name,
    start,
    GROUP_CONCAT(CASE WHEN cnt = 1 THEN 'S'
      WHEN cnt = 2 THEN 'D'
      ELSE STRING(cnt) END) AS SINGLETON_DOUBLETON,
    reference_bases,
    alternate_bases,
    GROUP_CONCAT(call.call_set_name) AS call.call_set_name,
    GROUP_CONCAT(genotype) AS genotype,
    SUM(num_samples_with_variant) AS num_samples_with_variant
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      alternate_bases,
      alt_num,
      call.call_set_name,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
      SUM(call.genotype == alt_num) WITHIN call AS cnt,
      COUNT(call.call_set_name) WITHIN RECORD AS num_samples_with_variant
    FROM (
        FLATTEN((
          SELECT
            reference_name,
            start,
            reference_bases,
            alternate_bases,
            POSITION(alternate_bases) AS alt_num,
            call.call_set_name,
            call.genotype,
          FROM
            [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
          # Optionally add a clause here to limit the query to a particular
          # region of the genome.
          #_WHERE_
          OMIT call IF EVERY(call.genotype = -1)
        ), alternate_bases)
        )
    OMIT RECORD IF alternate_bases IS NULL
    HAVING
      cnt > 0
      )
    GROUP EACH BY
      reference_name,
      start,
      reference_bases,
      alternate_bases
  HAVING
    num_samples_with_variant = 1
    )
GROUP BY
  call.call_set_name
ORDER BY
  private_variant_count DESC

Running query:   RUNNING  2.1s
Running query:   RUNNING  2.8s
Running query:   RUNNING  3.4s
Running query:   RUNNING  4.7s
Running query:   RUNNING  5.3s
Running query:   RUNNING  5.9s
Running query:   RUNNING  6.6s
```
Number of rows returned by this query: 5.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:43:20 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 554759 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 542170 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 534536 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 508429 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 501088 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result) +
  geom_point(aes(x=call_call_set_name, y=private_variant_count)) +
  theme(axis.text.x=if(nrow(result) <= 20)
    {element_text(angle = 90, hjust = 1)} else {element_blank()}) +
  xlab("Sample") +
  ylab("Number of Singletons") +
  ggtitle("Count of Singletons Per Genome")
```

<img src="figure/singletons-1.png" title="plot of chunk singletons" alt="plot of chunk singletons" style="display: block; margin: auto;" />

## Homozygosity Rate and Inbreeding Coefficient

For each genome, compare the expected and observed rates of homozygosity.


```r
result <- DisplayAndDispatchQuery("./sql/homozygous-variants.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Compute the expected and observed homozygosity rate for each individual.
SELECT
  call.call_set_name,
  O_HOM,
  ROUND(E_HOM, 2) as E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call.call_set_name,
    SUM(first_allele = second_allele) AS O_HOM,
    SUM(1.0 - (2.0 * freq * (1.0 - freq) * (called_allele_count / (called_allele_count - 1.0)))) AS E_HOM,
    COUNT(call.call_set_name) AS N_SITES,
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF((SUM(1 = call.genotype) > 0),
        SUM(call.genotype = 1)/SUM(call.genotype >= 0),
        -1)  WITHIN RECORD AS freq
    FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR (2 > COUNT(call.genotype))
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    call.call_set_name
    )
ORDER BY
  call.call_set_name

Running query:   RUNNING  2.7s
Running query:   RUNNING  3.4s
```
Number of rows returned by this query: 5.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:43:28 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 4412463 </td> <td align="right"> 5069576.52 </td> <td align="right"> 6580000 </td> <td align="right"> -0.44 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 4392441 </td> <td align="right"> 5038755.21 </td> <td align="right"> 6559521 </td> <td align="right"> -0.42 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 4350169 </td> <td align="right"> 5099971.91 </td> <td align="right"> 6600877 </td> <td align="right"> -0.50 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 4426399 </td> <td align="right"> 5065745.63 </td> <td align="right"> 6546838 </td> <td align="right"> -0.43 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 4385816 </td> <td align="right"> 5027331.32 </td> <td align="right"> 6537076 </td> <td align="right"> -0.42 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result) +
  geom_text(aes(x=O_HOM, y=E_HOM, label=call_call_set_name), hjust=-1, vjust=0) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="figure/homozygosity-1.png" title="plot of chunk homozygosity" alt="plot of chunk homozygosity" style="display: block; margin: auto;" />

## Sex Inference

For each genome, compare the gender from the sample information to the heterozygosity rate on the chromosome X calls.

```r
result <- DisplayAndDispatchQuery("./sql/gender-check.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Compute the the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the gender phenotype value is
# correct for each individual.
SELECT
  call.call_set_name,
  ROUND((het_RA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_het_alt_in_snvs,
  ROUND((hom_AA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_hom_alt_in_snvs,
  (hom_AA_count + het_RA_count + hom_RR_count) AS all_callable_sites,
  hom_AA_count,
  het_RA_count,
  hom_RR_count,
  (hom_AA_count + het_RA_count) AS all_snvs,
FROM
  (
  SELECT
    call.call_set_name,
    SUM(0 = first_allele
      AND 0 = second_allele) AS hom_RR_count,
    SUM(first_allele = second_allele AND first_allele > 0) AS hom_AA_count,
    SUM((first_allele != second_allele OR second_allele IS NULL)
      AND (first_allele > 0 OR second_allele > 0)) AS het_RA_count
  FROM (
    SELECT
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
    FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls]
    WHERE
      reference_name = 'chrX'
      AND start NOT BETWEEN 59999 AND 2699519
      AND start NOT BETWEEN 154931042 AND 155260559
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    call.call_set_name)
ORDER BY
  call.call_set_name
```
Number of rows returned by this query: 5.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:43:33 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> all_callable_sites </th> <th> hom_AA_count </th> <th> het_RA_count </th> <th> hom_RR_count </th> <th> all_snvs </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 166048 </td> <td align="right"> 71719 </td> <td align="right"> 2090 </td> <td align="right"> 92239 </td> <td align="right"> 73809 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 165523 </td> <td align="right"> 72908 </td> <td align="right"> 1950 </td> <td align="right"> 90665 </td> <td align="right"> 74858 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 0.59 </td> <td align="right"> 0.41 </td> <td align="right"> 167209 </td> <td align="right"> 44798 </td> <td align="right"> 64482 </td> <td align="right"> 57929 </td> <td align="right"> 109280 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 165883 </td> <td align="right"> 74881 </td> <td align="right"> 2165 </td> <td align="right"> 88837 </td> <td align="right"> 77046 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 165372 </td> <td align="right"> 73620 </td> <td align="right"> 2592 </td> <td align="right"> 89160 </td> <td align="right"> 76212 </td> </tr>
   </table>

Let's join this with the sample information:

```r
joinedResult <- inner_join(result, sampleInfo)
```

And visualize the results:

```r
ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender)) +
  theme(axis.text.x=if(nrow(result) <= 20)
    {element_text(angle = 90, hjust = 1)} else {element_blank()}) +
  xlab("Sample") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Heterozygosity Rate on the X Chromosome")
```

<img src="figure/gender-1.png" title="plot of chunk gender" alt="plot of chunk gender" style="display: block; margin: auto;" />

## Ethnicity Inference

For each genome, compare the ethncity from the sample information to the clustering in this analysis.

For this check, we:
* use the intersection of common variants found in both 1,000 Genomes phase 1 variants and Platinum Genomes
* compute PCA on those variants in common between the two data
* examine whether the individuals in Platinum Genomes cluster with other samples of the same ethnicity

Note that this `n^2` analysis is a cluster compute job instead of a BigQuery query.

This is a work-in-progress.  See https://github.com/elmer-garduno/spark-examples/tree/multiple_dataset_pca for the current state.

## Genome Similarity

Perform a simplistic similarity check on each pair of genomes to identify any mislabled or cross-contaminated samples.

Note that this `n^2` analysis is a cluster compute job instead of a BigQuery query.

### Results

This is hard-coded to Identity-By-State results for Platinum Genomes.

```r
ibs <- read.table("./data/platinum-genomes-ibs.tsv",
                  col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))
```

```
Warning in file(file, "rt"): cannot open file
'./data/platinum-genomes-ibs.tsv': No such file or directory
```

```
Error in file(file, "rt"): cannot open the connection
```

```r
ggplot(ibs) +
  geom_tile(aes(x=sample1, y=sample2, fill=ibsScore), colour="white") +
  scale_fill_gradient(low="white", high="steelblue",
                      na.value="black", trans="log",
                      guide=guide_colourbar(title= "IBS Score")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample 1") +
  ylab("Sample 2") +
  ggtitle("Identity By State (IBS) Heat Map")
```

```
Error in ggplot(ibs): object 'ibs' not found
```

### To Run the Cluster Compute Job

If you wish to run the Dataflow job, see the [dataflow-java README](https://github.com/googlegenomics/dataflow-java) for instructions to compile and run the job.
```
java -cp target/google-genomics-dataflow-v1beta2-0.2-SNAPSHOT.jar \
com.google.cloud.genomics.dataflow.pipelines.IdentityByState \
--project=YOUR-PROJECT \
--stagingLocation=gs://YOUR-BUCKET/staging \
--output=gs://YOUR-BUCKET/output/platinum-genomes-ibs.tsv \
--genomicsSecretsFile=/PATH/TO/YOUR/client_secrets.json \
--runner=DataflowPipelineRunner \
--numWorkers=40 \
--basesPerShard=1000000 \
--datasetId=3049512673186936334 \
--nonVariantSegments \
--allReferences
```

Note that there are several IBS calculators from which to choose.  Use the `--callSimilarityCalculatorFactory` to switch between them.

To run the job on a different dataset, change the variant set id for the `--datasetId` id parameter.  (Also, remove the `--nonVariantSegments` parameter if the data does not contain them.)

To gather the results into a single file:
```
gsutil cat gs://YOUR-BUCKET/output/platinum-genomes-ibs.tsv* | sort > platinum-genomes-ibs.tsv
```

# Removing Genomes from the Cohort

To remove a genome from a variant set in the Genomics API:
* See the [callsets delete](https://cloud.google.com/genomics/v1beta2/reference/callsets/delete) method.
* To delete a callset using a command line tool, see the the `deletecallset` command in [api-client-java](http://github.com/googlegenomics/api-client-java).

To only remove a genome from BigQuery only:
* Re-export the table to BigQuery using the `--call_set_id` flag on the `exportvariants` command in [api-client-java](http://github.com/googlegenomics/api-client-java) to list which callsets to _include_ in the export.

--------------------------------------------------------
_Next_: [Part 4: Variant-Level QC](./Variant-Level-QC.md)
