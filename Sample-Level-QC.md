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

*This code is made in collaboration with [Google Genomics](https://github.com/googlegenomics)*

# Part 3: Sample-Level QC





In Part 3 of the codelab, we perform some quality control analyses that could help to identify any problematic genomes that should be removed from the cohort before proceeding with further analysis.  The appropriate cut off thresholds will depend upon the input dataset and/or other factors.

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
tableReplacement <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs_no_tissue",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java2",
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data")
sampleData <- read.csv("http://storage.googleapis.com/genomics-public-data/platinum-genomes/other/platinum_genomes_sample_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

ibs <- read.table("./data/all-genomes-ibs.tsv",
                  col.names=c("sample1", "sample2", "ibsScore", "similar", "observed"))

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpCGIOnlyDataset.R")
```


* [Missingness Rate](#missingness-rate)
* [Singleton Rate](#singleton-rate)
* [Heterozygosity Rate and Inbreeding Coefficient](#homozygosity-rate-and-inbreeding-coefficient)
* [Sex Inference](#sex-inference)
* [Genotyping Concordance](#genotyping-concordance)
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
result <- DisplayAndDispatchQuery("./sql/missingness-sample-level.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
# Determine the level of missingness for each sample when compared to the hg19 reference genome
SELECT 
  g.sample_id AS sample_id,
  ROUND(((hg19.count - g.all_calls_count)/hg19.count), 3) AS missingness
FROM (
  SELECT
    call.call_set_name AS sample_id,
    ref_count + alt_count AS all_calls_count,
  FROM (
    SELECT
      call.call_set_name,
      SUM(IF(genotype = '0,0',
        (end - start),
        0)) AS ref_count,
      SUM(IF(genotype NOT IN ('0,0', '-1,-1'),
        1,
        0)) AS alt_count,
    FROM (
    SELECT
      call.call_set_name AS sample_id,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
    FROM
      [va_aaa_pilot_data.all_genomes_gvcfs_no_tissue])
    OMIT call IF call.call_set_name = 'LP6005243-DNA_A08'
    GROUP BY
      call.call_set_name)) AS g
  CROSS JOIN (
    SELECT 
      COUNT(Chr) AS count
    FROM 
      [google.com:biggene:test.hg19]) AS hg19
ORDER BY
  missingness DESC
Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.9s
Running query:   RUNNING  4.5s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.8s
```
Number of rows returned by this query: 458.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:05:24 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> missingness </th>  </tr>
  <tr> <td> LP6005038-DNA_C02 </td> <td align="right"> 0.06 </td> </tr>
  <tr> <td> LP6005692-DNA_D09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_D07 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_D02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_G10 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_G01 </td> <td align="right"> 0.05 </td> </tr>
   </table>

And visualizing the results:

```r
ggplot(result) +
  geom_point(aes(x=sample_id, y=missingness)) +
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
            [va_aaa_pilot_data.all_genomes_gvcfs_no_tissue]
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
```
Number of rows returned by this query: 459.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:05:28 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> LP6005051-DNA_D09 </td> <td align="right"> 500391 </td> </tr>
  <tr> <td> LP6005692-DNA_D05 </td> <td align="right"> 405391 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td align="right"> 171242 </td> </tr>
  <tr> <td> LP6005692-DNA_H01 </td> <td align="right"> 137990 </td> </tr>
  <tr> <td> LP6005144-DNA_D03 </td> <td align="right"> 136569 </td> </tr>
  <tr> <td> LP6005243-DNA_G07 </td> <td align="right"> 82186 </td> </tr>
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
      [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
    # Optionally add a clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    # Skip no calls and haploid sites
    OMIT call IF SOME(call.genotype < 0) OR (2 > COUNT(call.genotype)) OR (call.call_set_name = 'LP6005243-DNA_A08')
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
```

```
Error: Cannot query the cross product of repeated fields alternate_bases and call.call_set_name. 

query invalidQuery. Cannot query the cross product of repeated fields alternate_bases and call.call_set_name. 
```
Number of rows returned by this query: 459.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:05:31 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> LP6005051-DNA_D09 </td> <td align="right"> 500391 </td> </tr>
  <tr> <td> LP6005692-DNA_D05 </td> <td align="right"> 405391 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td align="right"> 171242 </td> </tr>
  <tr> <td> LP6005692-DNA_H01 </td> <td align="right"> 137990 </td> </tr>
  <tr> <td> LP6005144-DNA_D03 </td> <td align="right"> 136569 </td> </tr>
  <tr> <td> LP6005243-DNA_G07 </td> <td align="right"> 82186 </td> </tr>
   </table>

And visualizing the results:

```r
limits <- c(min(result$O_HOM, result$E_HOM),
            max(result$O_HOM, result$E_HOM))
```

```
Warning in min(result$O_HOM, result$E_HOM): no non-missing arguments to
min; returning Inf
```

```
Warning in max(result$O_HOM, result$E_HOM): no non-missing arguments to
max; returning -Inf
```

```r
ggplot(result) +
  geom_point(aes(x=O_HOM, y=E_HOM, label=call_call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits) + 
  scale_y_continuous(limits=limits) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

```
Error in eval(expr, envir, enclos): object 'O_HOM' not found
```

And with labels:

```r
ggplot(result) +
  geom_text(aes(x=O_HOM, y=E_HOM, label=call_call_set_name, hjust=0, vjust=0), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits, expand=c(0.05, 5)) +
  scale_y_continuous(limits=limits) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

```
Error in eval(expr, envir, enclos): object 'O_HOM' not found
```


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
      [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
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
Number of rows returned by this query: 478.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:05:34 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> all_callable_sites </th> <th> hom_AA_count </th> <th> het_RA_count </th> <th> hom_RR_count </th> <th> all_snvs </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 879381 </td> <td align="right"> 71719 </td> <td align="right"> 2090 </td> <td align="right"> 805572 </td> <td align="right"> 73809 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 873063 </td> <td align="right"> 72907 </td> <td align="right"> 1949 </td> <td align="right"> 798207 </td> <td align="right"> 74856 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 0.59 </td> <td align="right"> 0.41 </td> <td align="right"> 877397 </td> <td align="right"> 44797 </td> <td align="right"> 64482 </td> <td align="right"> 768118 </td> <td align="right"> 109279 </td> </tr>
  <tr> <td> LP6005038-DNA_A04 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 875855 </td> <td align="right"> 71427 </td> <td align="right"> 2215 </td> <td align="right"> 802213 </td> <td align="right"> 73642 </td> </tr>
  <tr> <td> LP6005038-DNA_A05 </td> <td align="right"> 0.58 </td> <td align="right"> 0.42 </td> <td align="right"> 875866 </td> <td align="right"> 44537 </td> <td align="right"> 61520 </td> <td align="right"> 769809 </td> <td align="right"> 106057 </td> </tr>
  <tr> <td> LP6005038-DNA_A06 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 874831 </td> <td align="right"> 76007 </td> <td align="right"> 2232 </td> <td align="right"> 796592 </td> <td align="right"> 78239 </td> </tr>
   </table>

Let's join this with the sample information:

```r
joinedResult <- inner_join(result, sampleInfo)
```

And visualize the results:

```r
ggplot(joinedResult) +
  geom_boxplot(aes(x=gender, y=perct_het_alt_in_snvs, fill=gender)) +
  scale_y_continuous() +
  xlab("Gender") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Box Plot: Heterozygosity Rate on the X Chromosome")
```

<img src="figure/gender-boxplot-1.png" title="plot of chunk gender-boxplot" alt="plot of chunk gender-boxplot" style="display: block; margin: auto;" />


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

## Genotyping Concordance

We next want to look at the concordance between SNPs called from the sequencing data and those called through the use genotyping.  This allows us to identify samples that may have been mixed up in the laboratory.


```r
concordanceResult <- DisplayAndDispatchQuery("./sql/genotyping-concordance.sql",
                                  project=project,
                                  replacements=tableReplacement)
```

```
SELECT
  sample_id,
  calls_in_common,
  identical_calls,
  (identical_calls/calls_in_common) AS concordance
FROM (
  SELECT 
    sample_id,
    COUNT(seq_genotype) AS calls_in_common,
    SUM(IF(seq_genotype = gen_genotype, 1, 0)) AS identical_calls,
  FROM (
    SELECT
      seq.sample_id AS sample_id,
      seq.reference_name AS reference_name,
      seq.start AS start,
      seq.end AS end,
      seq.genotype AS seq_genotype,
      gen.genotype AS gen_genotype,
    FROM (
      SELECT
        sample_id,
        reference_name,
        start,
        end,
        genotype,
        bin,
      FROM js(
        (SELECT
          call.call_set_name,
          reference_name,
          start,
          end,
          call.genotype,
          reference_bases,
          GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alts,
          COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        FROM
          [va_aaa_pilot_data.all_genomes_gvcfs_no_tissue]
         #_WHERE_
        OMIT 
          call IF EVERY (call.genotype < 0)
        HAVING 
          num_alts <= 1
          AND reference_bases IN ('A','C','G','T')
          AND (alts IS null
            OR LENGTH(alts) <= 1)
        ),
        // Start javascript function
        // Input Columns
        call.call_set_name, reference_name, start, end, call.genotype,
        // Output Schema
        "[{name: 'sample_id', type: 'string'},
        {name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'end', type: 'integer'},
        {name: 'genotype', type: 'string'},
        {name: 'bin', type: 'integer'}]",
        // Function
        "function(r, emit) {
          for (c of r.call) {
            var binSize = 5000;
            var startBin = Math.floor(r.start/binSize);
            var endBin = Math.floor(r.end/binSize);
            var genotype = JSON.stringify(c.genotype.sort());
            for (var bin = startBin; bin <= endBin; bin++){
              emit({
                sample_id: c.call_set_name,
                reference_name: r.reference_name,
                start: r.start,
                end: r.end,
                genotype: genotype,
                bin: bin,
              })
            }
          }
        }")) AS seq
JOIN EACH (
  SELECT
    sample_id,
    reference_name,
    start,
    end,
    genotype,
    bin,
  FROM js(
    (SELECT
      call.call_set_name,
      reference_name,
      start,
      end,
      call.genotype,
    FROM
      [va_aaa_pilot_data.genotyping_data]
      OMIT call IF EVERY (call.genotype < 0)       
    ),
    // Start javascript function
    // Input Columns
    call.call_set_name, reference_name, start, end, call.genotype,
    // Output Schema
    "[{name: 'sample_id', type: 'string'},
    {name: 'reference_name', type: 'string'},
    {name: 'start', type: 'integer'},
    {name: 'end', type: 'integer'},
    {name: 'genotype', type: 'string'},
    {name: 'bin', type: 'integer'}]",
    // Function
    "function(r, emit) {
      for (c of r.call) {
        var binSize = 5000;
        var bin = Math.floor(r.start/binSize);
        var genotype = JSON.stringify(c.genotype.sort());
        emit({
          sample_id: c.call_set_name,
          reference_name: r.reference_name,
          start: r.start,
          end: r.end,
          genotype: genotype,
          bin: bin,
        })
      }
    }")) AS gen
ON
  seq.sample_id = gen.sample_id
  AND seq.reference_name = gen.reference_name
  AND seq.bin = gen.bin
WHERE
  seq.start <= gen.start
  AND seq.end >= gen.end )
GROUP BY 
  sample_id)
Running query:   RUNNING  2.2s
Running query:   RUNNING  2.8s
Running query:   RUNNING  3.5s
Running query:   RUNNING  4.1s
Running query:   RUNNING  4.8s
Running query:   RUNNING  5.4s
Running query:   RUNNING  6.0s
Running query:   RUNNING  6.7s
Running query:   RUNNING  7.3s
Running query:   RUNNING  7.9s
Running query:   RUNNING  8.6s
Running query:   RUNNING  9.2s
Running query:   RUNNING  9.9s
Running query:   RUNNING 10.5s
Running query:   RUNNING 11.1s
Running query:   RUNNING 11.7s
Running query:   RUNNING 12.4s
Running query:   RUNNING 13.1s
Running query:   RUNNING 13.7s
Running query:   RUNNING 14.4s
Running query:   RUNNING 15.0s
Running query:   RUNNING 15.6s
Running query:   RUNNING 16.3s
Running query:   RUNNING 16.9s
Running query:   RUNNING 17.5s
Running query:   RUNNING 18.1s
Running query:   RUNNING 18.8s
Running query:   RUNNING 19.4s
Running query:   RUNNING 20.1s
Running query:   RUNNING 20.7s
Running query:   RUNNING 21.3s
Running query:   RUNNING 21.9s
Running query:   RUNNING 22.6s
Running query:   RUNNING 23.2s
Running query:   RUNNING 23.8s
Running query:   RUNNING 24.5s
Running query:   RUNNING 25.1s
Running query:   RUNNING 25.7s
Running query:   RUNNING 26.4s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.7s
Running query:   RUNNING 28.3s
Running query:   RUNNING 29.0s
Running query:   RUNNING 29.6s
Running query:   RUNNING 30.3s
Running query:   RUNNING 30.9s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.1s
Running query:   RUNNING 32.8s
Running query:   RUNNING 33.4s
Running query:   RUNNING 34.0s
Running query:   RUNNING 34.6s
Running query:   RUNNING 35.3s
Running query:   RUNNING 35.9s
Running query:   RUNNING 36.5s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.8s
Running query:   RUNNING 38.4s
Running query:   RUNNING 39.0s
Running query:   RUNNING 39.7s
Running query:   RUNNING 40.3s
Running query:   RUNNING 41.0s
Running query:   RUNNING 41.7s
Running query:   RUNNING 42.3s
Running query:   RUNNING 43.0s
Running query:   RUNNING 43.6s
Running query:   RUNNING 44.2s
Running query:   RUNNING 44.9s
Running query:   RUNNING 45.5s
Running query:   RUNNING 46.1s
Running query:   RUNNING 46.8s
Running query:   RUNNING 47.4s
Running query:   RUNNING 48.0s
Running query:   RUNNING 48.7s
Running query:   RUNNING 49.3s
Running query:   RUNNING 49.9s
Running query:   RUNNING 50.6s
Running query:   RUNNING 51.2s
Running query:   RUNNING 51.8s
Running query:   RUNNING 52.5s
Running query:   RUNNING 53.1s
Running query:   RUNNING 54.0s
Running query:   RUNNING 54.6s
Running query:   RUNNING 55.3s
Running query:   RUNNING 55.9s
Running query:   RUNNING 56.5s
Running query:   RUNNING 57.2s
Running query:   RUNNING 57.9s
Running query:   RUNNING 58.5s
Running query:   RUNNING 59.2s
Running query:   RUNNING 60.1s
Running query:   RUNNING 60.7s
Running query:   RUNNING 61.3s
Running query:   RUNNING 62.0s
Running query:   RUNNING 62.6s
Running query:   RUNNING 63.2s
Running query:   RUNNING 63.8s
Running query:   RUNNING 64.5s
Running query:   RUNNING 65.1s
Running query:   RUNNING 65.8s
Running query:   RUNNING 66.4s
Running query:   RUNNING 67.1s
Running query:   RUNNING 67.7s
Running query:   RUNNING 68.3s
Running query:   RUNNING 69.0s
Running query:   RUNNING 69.7s
Running query:   RUNNING 70.3s
Running query:   RUNNING 71.0s
Running query:   RUNNING 71.7s
Running query:   RUNNING 72.3s
Running query:   RUNNING 72.9s
Running query:   RUNNING 73.6s
Running query:   RUNNING 74.3s
Running query:   RUNNING 74.9s
Running query:   RUNNING 75.5s
Running query:   RUNNING 76.2s
Running query:   RUNNING 76.8s
Running query:   RUNNING 77.5s
Running query:   RUNNING 78.1s
Running query:   RUNNING 78.7s
Running query:   RUNNING 79.4s
Running query:   RUNNING 80.0s
Running query:   RUNNING 80.6s
Running query:   RUNNING 81.3s
Running query:   RUNNING 81.9s
Running query:   RUNNING 82.5s
Running query:   RUNNING 83.2s
Running query:   RUNNING 83.8s
Running query:   RUNNING 84.4s
Running query:   RUNNING 85.1s
Running query:   RUNNING 85.7s
Running query:   RUNNING 86.3s
Running query:   RUNNING 87.0s
Running query:   RUNNING 87.7s
Running query:   RUNNING 88.3s
Running query:   RUNNING 88.9s
Running query:   RUNNING 89.6s
Running query:   RUNNING 90.2s
Running query:   RUNNING 90.8s
Running query:   RUNNING 91.5s
Running query:   RUNNING 92.1s
Running query:   RUNNING 92.8s
Running query:   RUNNING 93.4s
Running query:   RUNNING 94.0s
Running query:   RUNNING 94.6s
Running query:   RUNNING 95.3s
Running query:   RUNNING 95.9s
Running query:   RUNNING 96.6s
Running query:   RUNNING 97.2s
Running query:   RUNNING 97.9s
Running query:   RUNNING 98.6s
Running query:   RUNNING 99.2s
Running query:   RUNNING 99.8s
Running query:   RUNNING 100.5s
Running query:   RUNNING 101.1s
Running query:   RUNNING 101.8s
Running query:   RUNNING 102.4s
Running query:   RUNNING 103.0s
Running query:   RUNNING 103.7s
Running query:   RUNNING 104.3s
Running query:   RUNNING 104.9s
Running query:   RUNNING 105.5s
Running query:   RUNNING 106.2s
Running query:   RUNNING 106.8s
Running query:   RUNNING 107.4s
Running query:   RUNNING 108.1s
Running query:   RUNNING 108.7s
Running query:   RUNNING 109.3s
Running query:   RUNNING 109.9s
Running query:   RUNNING 110.6s
Running query:   RUNNING 111.2s
Running query:   RUNNING 111.8s
Running query:   RUNNING 112.5s
Running query:   RUNNING 113.1s
Running query:   RUNNING 113.8s
Running query:   RUNNING 114.4s
Running query:   RUNNING 115.0s
Running query:   RUNNING 115.7s
Running query:   RUNNING 116.3s
Running query:   RUNNING 117.0s
Running query:   RUNNING 117.6s
Running query:   RUNNING 118.3s
Running query:   RUNNING 118.9s
Running query:   RUNNING 119.5s
Running query:   RUNNING 120.2s
Running query:   RUNNING 120.8s
Running query:   RUNNING 121.4s
Running query:   RUNNING 122.0s
Running query:   RUNNING 122.7s
Running query:   RUNNING 123.3s
Running query:   RUNNING 123.9s
Running query:   RUNNING 124.6s
Running query:   RUNNING 125.2s
Running query:   RUNNING 125.9s
Running query:   RUNNING 126.5s
Running query:   RUNNING 127.2s
Running query:   RUNNING 127.8s
Running query:   RUNNING 128.5s
Running query:   RUNNING 129.1s
Running query:   RUNNING 129.7s
Running query:   RUNNING 130.4s
Running query:   RUNNING 131.0s
Running query:   RUNNING 131.6s
Running query:   RUNNING 132.3s
Running query:   RUNNING 132.9s
Running query:   RUNNING 133.5s
Running query:   RUNNING 134.4s
Running query:   RUNNING 135.0s
Running query:   RUNNING 135.6s
Running query:   RUNNING 136.2s
Running query:   RUNNING 136.9s
Running query:   RUNNING 137.5s
Running query:   RUNNING 138.2s
Running query:   RUNNING 138.8s
Running query:   RUNNING 139.4s
Running query:   RUNNING 140.1s
Running query:   RUNNING 140.7s
Running query:   RUNNING 141.3s
Running query:   RUNNING 142.0s
Running query:   RUNNING 142.6s
Running query:   RUNNING 143.2s
Running query:   RUNNING 143.9s
Running query:   RUNNING 144.5s
Running query:   RUNNING 145.1s
Running query:   RUNNING 145.8s
Running query:   RUNNING 146.4s
Running query:   RUNNING 147.1s
Running query:   RUNNING 147.7s
Running query:   RUNNING 148.3s
Running query:   RUNNING 149.0s
Running query:   RUNNING 149.6s
Running query:   RUNNING 150.2s
Running query:   RUNNING 150.9s
Running query:   RUNNING 151.5s
Running query:   RUNNING 152.2s
Running query:   RUNNING 152.8s
Running query:   RUNNING 153.4s
Running query:   RUNNING 154.1s
Running query:   RUNNING 154.7s
Running query:   RUNNING 155.4s
Running query:   RUNNING 156.0s
Running query:   RUNNING 156.6s
Running query:   RUNNING 157.3s
Running query:   RUNNING 157.9s
Running query:   RUNNING 158.5s
Running query:   RUNNING 159.2s
Running query:   RUNNING 159.8s
Running query:   RUNNING 160.4s
Running query:   RUNNING 161.1s
Running query:   RUNNING 161.7s
Running query:   RUNNING 162.3s
Running query:   RUNNING 163.0s
Running query:   RUNNING 163.6s
Running query:   RUNNING 164.3s
Running query:   RUNNING 165.2s
Running query:   RUNNING 165.8s
Running query:   RUNNING 166.5s
Running query:   RUNNING 167.1s
Running query:   RUNNING 167.7s
Running query:   RUNNING 168.4s
Running query:   RUNNING 169.1s
Running query:   RUNNING 169.7s
Running query:   RUNNING 170.4s
Running query:   RUNNING 171.0s
Running query:   RUNNING 171.6s
Running query:   RUNNING 172.2s
Running query:   RUNNING 172.9s
Running query:   RUNNING 173.6s
Running query:   RUNNING 174.2s
Running query:   RUNNING 175.6s
Running query:   RUNNING 176.2s
Running query:   RUNNING 176.8s
Running query:   RUNNING 177.5s
Running query:   RUNNING 178.1s
Running query:   RUNNING 178.7s
Running query:   RUNNING 179.4s
Running query:   RUNNING 180.1s
Running query:   RUNNING 180.7s
Running query:   RUNNING 181.3s
Running query:   RUNNING 182.0s
Running query:   RUNNING 182.7s
Running query:   RUNNING 183.3s
Running query:   RUNNING 183.9s
Running query:   RUNNING 184.6s
Running query:   RUNNING 185.2s
Running query:   RUNNING 185.8s
Running query:   RUNNING 186.5s
Running query:   RUNNING 187.1s
Running query:   RUNNING 187.8s
Running query:   RUNNING 188.4s
Running query:   RUNNING 189.0s
Running query:   RUNNING 189.6s
Running query:   RUNNING 190.3s
Running query:   RUNNING 190.9s
Running query:   RUNNING 191.5s
Running query:   RUNNING 192.2s
Running query:   RUNNING 192.8s
Running query:   RUNNING 193.4s
Running query:   RUNNING 194.0s
Running query:   RUNNING 194.7s
Running query:   RUNNING 195.3s
Running query:   RUNNING 196.0s
Running query:   RUNNING 196.6s
Running query:   RUNNING 197.2s
Running query:   RUNNING 197.9s
Running query:   RUNNING 198.5s
Running query:   RUNNING 199.1s
Running query:   RUNNING 199.8s
Running query:   RUNNING 200.4s
Running query:   RUNNING 201.0s
Running query:   RUNNING 201.7s
Running query:   RUNNING 202.3s
Running query:   RUNNING 202.9s
Running query:   RUNNING 203.6s
Running query:   RUNNING 204.2s
Running query:   RUNNING 204.8s
Running query:   RUNNING 205.5s
Running query:   RUNNING 206.1s
Running query:   RUNNING 206.8s
Running query:   RUNNING 207.4s
Running query:   RUNNING 208.0s
Running query:   RUNNING 208.6s
Running query:   RUNNING 209.3s
Running query:   RUNNING 209.9s
Running query:   RUNNING 210.6s
Running query:   RUNNING 211.3s
Running query:   RUNNING 211.9s
Running query:   RUNNING 212.5s
Running query:   RUNNING 213.2s
Running query:   RUNNING 213.9s
Running query:   RUNNING 214.5s
Running query:   RUNNING 215.2s
Running query:   RUNNING 215.8s
Running query:   RUNNING 216.4s
Running query:   RUNNING 217.0s
Running query:   RUNNING 217.7s
Running query:   RUNNING 218.3s
Running query:   RUNNING 219.0s
Running query:   RUNNING 219.7s
Running query:   RUNNING 220.3s
Running query:   RUNNING 220.9s
Running query:   RUNNING 221.6s
Running query:   RUNNING 222.2s
Running query:   RUNNING 222.9s
Running query:   RUNNING 223.5s
Running query:   RUNNING 224.2s
Running query:   RUNNING 224.8s
Running query:   RUNNING 225.4s
Running query:   RUNNING 226.0s
Running query:   RUNNING 226.7s
Running query:   RUNNING 227.3s
Running query:   RUNNING 228.0s
Running query:   RUNNING 228.6s
Running query:   RUNNING 229.2s
Running query:   RUNNING 229.9s
Running query:   RUNNING 230.5s
Running query:   RUNNING 231.1s
Running query:   RUNNING 231.8s
Running query:   RUNNING 232.4s
Running query:   RUNNING 233.1s
Running query:   RUNNING 233.7s
Running query:   RUNNING 234.3s
Running query:   RUNNING 235.0s
Running query:   RUNNING 235.6s
Running query:   RUNNING 236.3s
Running query:   RUNNING 236.9s
Running query:   RUNNING 237.6s
Running query:   RUNNING 238.2s
Running query:   RUNNING 238.8s
Running query:   RUNNING 239.5s
Running query:   RUNNING 240.1s
Running query:   RUNNING 240.8s
Running query:   RUNNING 241.5s
Running query:   RUNNING 242.1s
Running query:   RUNNING 242.7s
Running query:   RUNNING 243.4s
Running query:   RUNNING 244.0s
Running query:   RUNNING 244.6s
Running query:   RUNNING 245.3s
Running query:   RUNNING 245.9s
Running query:   RUNNING 246.5s
Running query:   RUNNING 247.2s
Running query:   RUNNING 247.8s
Running query:   RUNNING 248.5s
Running query:   RUNNING 249.1s
Running query:   RUNNING 249.7s
Running query:   RUNNING 250.4s
Running query:   RUNNING 251.0s
Running query:   RUNNING 251.7s
Running query:   RUNNING 252.3s
Running query:   RUNNING 253.0s
Running query:   RUNNING 253.6s
Running query:   RUNNING 254.2s
Running query:   RUNNING 254.8s
Running query:   RUNNING 255.5s
Running query:   RUNNING 256.1s
Running query:   RUNNING 256.8s
Running query:   RUNNING 257.5s
Running query:   RUNNING 258.1s
Running query:   RUNNING 258.7s
Running query:   RUNNING 259.4s
Running query:   RUNNING 260.0s
Running query:   RUNNING 260.6s
Running query:   RUNNING 261.3s
Running query:   RUNNING 261.9s
Running query:   RUNNING 262.6s
Running query:   RUNNING 263.2s
Running query:   RUNNING 263.9s
Running query:   RUNNING 264.5s
Running query:   RUNNING 265.2s
Running query:   RUNNING 265.8s
Running query:   RUNNING 266.5s
Running query:   RUNNING 267.1s
Running query:   RUNNING 267.8s
Running query:   RUNNING 268.4s
Running query:   RUNNING 269.1s
Running query:   RUNNING 269.7s
Running query:   RUNNING 270.3s
Running query:   RUNNING 271.0s
Running query:   RUNNING 271.6s
Running query:   RUNNING 272.2s
Running query:   RUNNING 272.9s
Running query:   RUNNING 273.5s
Running query:   RUNNING 274.1s
Running query:   RUNNING 274.7s
Running query:   RUNNING 275.4s
Running query:   RUNNING 276.0s
Running query:   RUNNING 276.7s
Running query:   RUNNING 277.3s
Running query:   RUNNING 278.0s
Running query:   RUNNING 278.6s
Running query:   RUNNING 279.2s
Running query:   RUNNING 279.9s
Running query:   RUNNING 280.5s
Running query:   RUNNING 281.2s
Running query:   RUNNING 281.8s
Running query:   RUNNING 282.5s
Running query:   RUNNING 283.2s
Running query:   RUNNING 283.8s
Running query:   RUNNING 284.5s
Running query:   RUNNING 285.1s
Running query:   RUNNING 285.8s
Running query:   RUNNING 286.5s
Running query:   RUNNING 287.1s
Running query:   RUNNING 287.8s
Running query:   RUNNING 288.4s
Running query:   RUNNING 289.1s
Running query:   RUNNING 289.8s
Running query:   RUNNING 290.5s
Running query:   RUNNING 291.1s
Running query:   RUNNING 291.8s
Running query:   RUNNING 292.5s
Running query:   RUNNING 293.1s
Running query:   RUNNING 293.8s
Running query:   RUNNING 294.4s
Running query:   RUNNING 295.1s
Running query:   RUNNING 295.8s
Running query:   RUNNING 296.5s
Running query:   RUNNING 297.1s
Running query:   RUNNING 297.8s
Running query:   RUNNING 298.4s
Running query:   RUNNING 299.1s
Running query:   RUNNING 299.8s
Running query:   RUNNING 300.5s
Running query:   RUNNING 301.1s
Running query:   RUNNING 301.8s
Running query:   RUNNING 302.5s
Running query:   RUNNING 303.1s
Running query:   RUNNING 303.8s
Running query:   RUNNING 304.4s
Running query:   RUNNING 305.1s
Running query:   RUNNING 305.8s
Running query:   RUNNING 306.5s
Running query:   RUNNING 307.1s
Running query:   RUNNING 307.8s
Running query:   RUNNING 308.5s
Running query:   RUNNING 309.2s
Running query:   RUNNING 309.8s
Running query:   RUNNING 310.5s
Running query:   RUNNING 311.2s
Running query:   RUNNING 311.9s
Running query:   RUNNING 312.5s
Running query:   RUNNING 313.1s
Running query:   RUNNING 313.8s
Running query:   RUNNING 314.4s
Running query:   RUNNING 315.1s
Running query:   RUNNING 315.8s
Running query:   RUNNING 316.4s
Running query:   RUNNING 317.1s
Running query:   RUNNING 317.7s
Running query:   RUNNING 318.4s
Running query:   RUNNING 319.0s
Running query:   RUNNING 319.8s
Running query:   RUNNING 320.4s
Running query:   RUNNING 321.0s
Running query:   RUNNING 321.7s
Running query:   RUNNING 322.3s
Running query:   RUNNING 323.0s
Running query:   RUNNING 323.6s
Running query:   RUNNING 324.3s
Running query:   RUNNING 324.9s
Running query:   RUNNING 325.6s
Running query:   RUNNING 326.2s
Running query:   RUNNING 326.9s
Running query:   RUNNING 327.5s
Running query:   RUNNING 328.1s
Running query:   RUNNING 328.8s
Running query:   RUNNING 329.4s
Running query:   RUNNING 330.1s
Running query:   RUNNING 330.7s
Running query:   RUNNING 331.3s
Running query:   RUNNING 332.0s
Running query:   RUNNING 332.6s
Running query:   RUNNING 333.3s
Running query:   RUNNING 333.9s
Running query:   RUNNING 334.6s
Running query:   RUNNING 335.2s
Running query:   RUNNING 335.9s
Running query:   RUNNING 336.5s
Running query:   RUNNING 337.1s
Running query:   RUNNING 337.8s
Running query:   RUNNING 338.4s
Running query:   RUNNING 339.1s
Running query:   RUNNING 339.7s
Running query:   RUNNING 340.4s
Running query:   RUNNING 341.0s
Running query:   RUNNING 341.7s
Running query:   RUNNING 342.3s
Running query:   RUNNING 342.9s
Running query:   RUNNING 343.6s
Running query:   RUNNING 344.3s
Running query:   RUNNING 344.9s
Running query:   RUNNING 345.5s
Running query:   RUNNING 346.2s
Running query:   RUNNING 346.8s
Running query:   RUNNING 347.5s
Running query:   RUNNING 348.1s
Running query:   RUNNING 348.7s
Running query:   RUNNING 349.4s
Running query:   RUNNING 350.0s
Running query:   RUNNING 350.6s
Running query:   RUNNING 351.3s
Running query:   RUNNING 351.9s
Running query:   RUNNING 352.5s
Running query:   RUNNING 353.1s
Running query:   RUNNING 353.8s
Running query:   RUNNING 354.4s
Running query:   RUNNING 355.0s
Running query:   RUNNING 355.7s
Running query:   RUNNING 356.3s
Running query:   RUNNING 357.0s
Running query:   RUNNING 357.6s
Running query:   RUNNING 358.3s
Running query:   RUNNING 358.9s
Running query:   RUNNING 359.5s
Running query:   RUNNING 360.2s
Running query:   RUNNING 360.9s
Running query:   RUNNING 361.5s
Running query:   RUNNING 362.1s
Running query:   RUNNING 362.7s
Running query:   RUNNING 363.4s
Running query:   RUNNING 364.1s
Running query:   RUNNING 364.7s
Running query:   RUNNING 365.3s
Running query:   RUNNING 366.0s
Running query:   RUNNING 366.6s
Running query:   RUNNING 367.3s
Running query:   RUNNING 367.9s
Running query:   RUNNING 368.6s
Running query:   RUNNING 369.2s
Running query:   RUNNING 369.8s
Running query:   RUNNING 370.5s
Running query:   RUNNING 371.1s
Running query:   RUNNING 371.8s
Running query:   RUNNING 372.4s
Running query:   RUNNING 373.1s
Running query:   RUNNING 373.7s
Running query:   RUNNING 374.4s
Running query:   RUNNING 375.0s
Running query:   RUNNING 375.6s
Running query:   RUNNING 376.3s
Running query:   RUNNING 376.9s
Running query:   RUNNING 377.6s
Running query:   RUNNING 378.2s
Running query:   RUNNING 378.8s
Running query:   RUNNING 379.5s
Running query:   RUNNING 380.2s
Running query:   RUNNING 380.9s
Running query:   RUNNING 381.5s
Running query:   RUNNING 382.2s
Running query:   RUNNING 382.8s
Running query:   RUNNING 383.4s
Running query:   RUNNING 384.1s
Running query:   RUNNING 384.7s
Running query:   RUNNING 385.4s
Running query:   RUNNING 386.0s
Running query:   RUNNING 386.6s
Running query:   RUNNING 387.3s
Running query:   RUNNING 387.9s
Running query:   RUNNING 388.6s
Running query:   RUNNING 389.2s
Running query:   RUNNING 389.9s
Running query:   RUNNING 390.5s
Running query:   RUNNING 391.2s
Running query:   RUNNING 391.8s
Running query:   RUNNING 392.4s
Running query:   RUNNING 393.1s
Running query:   RUNNING 393.7s
Running query:   RUNNING 394.4s
Running query:   RUNNING 395.0s
Running query:   RUNNING 395.6s
Running query:   RUNNING 396.3s
Running query:   RUNNING 396.9s
Running query:   RUNNING 397.5s
Running query:   RUNNING 398.2s
Running query:   RUNNING 398.8s
Running query:   RUNNING 399.4s
Running query:   RUNNING 400.1s
Running query:   RUNNING 400.8s
Running query:   RUNNING 401.4s
Running query:   RUNNING 402.0s
Running query:   RUNNING 402.7s
Running query:   RUNNING 403.3s
Running query:   RUNNING 403.9s
Running query:   RUNNING 404.6s
Running query:   RUNNING 405.2s
Running query:   RUNNING 405.9s
Running query:   RUNNING 406.5s
Running query:   RUNNING 407.2s
Running query:   RUNNING 407.8s
Running query:   RUNNING 408.4s
Running query:   RUNNING 409.1s
Running query:   RUNNING 409.7s
Running query:   RUNNING 410.4s
Running query:   RUNNING 411.0s
Running query:   RUNNING 411.7s
Running query:   RUNNING 412.3s
Running query:   RUNNING 412.9s
Running query:   RUNNING 413.6s
Running query:   RUNNING 414.2s
Running query:   RUNNING 414.9s
Running query:   RUNNING 415.5s
Running query:   RUNNING 416.1s
Running query:   RUNNING 416.8s
Running query:   RUNNING 417.4s
Running query:   RUNNING 418.0s
Running query:   RUNNING 418.7s
Running query:   RUNNING 419.3s
Running query:   RUNNING 420.0s
Running query:   RUNNING 420.6s
Running query:   RUNNING 421.2s
Running query:   RUNNING 421.9s
Running query:   RUNNING 422.5s
Running query:   RUNNING 423.1s
Running query:   RUNNING 423.8s
Running query:   RUNNING 424.5s
Running query:   RUNNING 425.2s
Running query:   RUNNING 425.8s
Running query:   RUNNING 426.5s
Running query:   RUNNING 427.2s
Running query:   RUNNING 427.8s
Running query:   RUNNING 428.5s
Running query:   RUNNING 429.1s
Running query:   RUNNING 429.8s
Running query:   RUNNING 430.4s
Running query:   RUNNING 431.1s
Running query:   RUNNING 431.7s
Running query:   RUNNING 432.4s
Running query:   RUNNING 433.0s
Running query:   RUNNING 433.6s
Running query:   RUNNING 434.3s
Running query:   RUNNING 434.9s
Running query:   RUNNING 435.6s
Running query:   RUNNING 436.2s
Running query:   RUNNING 436.8s
Running query:   RUNNING 437.5s
Running query:   RUNNING 438.1s
Running query:   RUNNING 438.8s
Running query:   RUNNING 439.4s
Running query:   RUNNING 440.1s
Running query:   RUNNING 440.7s
Running query:   RUNNING 441.4s
Running query:   RUNNING 442.0s
Running query:   RUNNING 442.6s
Running query:   RUNNING 443.3s
Running query:   RUNNING 443.9s
Running query:   RUNNING 444.6s
Running query:   RUNNING 445.3s
Running query:   RUNNING 445.9s
Running query:   RUNNING 446.5s
Running query:   RUNNING 447.1s
Running query:   RUNNING 447.8s
Running query:   RUNNING 448.4s
Running query:   RUNNING 449.1s
Running query:   RUNNING 449.7s
Running query:   RUNNING 450.3s
Running query:   RUNNING 451.0s
Running query:   RUNNING 451.6s
Running query:   RUNNING 452.2s
Running query:   RUNNING 452.9s
Running query:   RUNNING 453.5s
Running query:   RUNNING 454.2s
Running query:   RUNNING 454.8s
Running query:   RUNNING 455.4s
Running query:   RUNNING 456.1s
Running query:   RUNNING 456.7s
Running query:   RUNNING 457.4s
Running query:   RUNNING 458.0s
Running query:   RUNNING 458.6s
Running query:   RUNNING 459.3s
Running query:   RUNNING 459.9s
Running query:   RUNNING 460.6s
Running query:   RUNNING 461.2s
Running query:   RUNNING 461.8s
Running query:   RUNNING 462.5s
Running query:   RUNNING 463.1s
Running query:   RUNNING 463.8s
Running query:   RUNNING 464.4s
Running query:   RUNNING 465.0s
Running query:   RUNNING 465.8s
Running query:   RUNNING 466.5s
Running query:   RUNNING 467.2s
Running query:   RUNNING 467.8s
Running query:   RUNNING 468.5s
Running query:   RUNNING 469.1s
Running query:   RUNNING 469.7s
Running query:   RUNNING 470.4s
Running query:   RUNNING 471.0s
Running query:   RUNNING 471.6s
Running query:   RUNNING 472.3s
Running query:   RUNNING 472.9s
Running query:   RUNNING 473.6s
Running query:   RUNNING 474.2s
Running query:   RUNNING 474.8s
Running query:   RUNNING 475.5s
Running query:   RUNNING 476.1s
Running query:   RUNNING 476.8s
Running query:   RUNNING 477.4s
Running query:   RUNNING 478.0s
Running query:   RUNNING 478.7s
Running query:   RUNNING 479.3s
Running query:   RUNNING 480.0s
Running query:   RUNNING 480.6s
Running query:   RUNNING 481.3s
Running query:   RUNNING 481.9s
Running query:   RUNNING 482.5s
Running query:   RUNNING 483.2s
Running query:   RUNNING 483.8s
Running query:   RUNNING 484.5s
Running query:   RUNNING 485.1s
Running query:   RUNNING 485.8s
Running query:   RUNNING 486.4s
Running query:   RUNNING 487.0s
Running query:   RUNNING 487.7s
Running query:   RUNNING 488.3s
Running query:   RUNNING 488.9s
Running query:   RUNNING 489.6s
Running query:   RUNNING 490.2s
Running query:   RUNNING 490.9s
Running query:   RUNNING 491.5s
Running query:   RUNNING 492.1s
Running query:   RUNNING 492.8s
Running query:   RUNNING 493.4s
Running query:   RUNNING 494.1s
Running query:   RUNNING 494.7s
Running query:   RUNNING 495.4s
Running query:   RUNNING 496.0s
Running query:   RUNNING 496.7s
Running query:   RUNNING 497.3s
Running query:   RUNNING 498.2s
Running query:   RUNNING 498.8s
Running query:   RUNNING 499.5s
Running query:   RUNNING 500.1s
Running query:   RUNNING 500.7s
Running query:   RUNNING 501.4s
Running query:   RUNNING 502.1s
Running query:   RUNNING 502.7s
Running query:   RUNNING 503.4s
Running query:   RUNNING 504.0s
Running query:   RUNNING 504.7s
Running query:   RUNNING 505.4s
Running query:   RUNNING 506.0s
Running query:   RUNNING 506.7s
Running query:   RUNNING 507.3s
Running query:   RUNNING 507.9s
Running query:   RUNNING 508.6s
Running query:   RUNNING 509.2s
Running query:   RUNNING 509.9s
Running query:   RUNNING 510.5s
Running query:   RUNNING 511.1s
Running query:   RUNNING 511.8s
Running query:   RUNNING 512.5s
Running query:   RUNNING 513.1s
Running query:   RUNNING 513.8s
Running query:   RUNNING 514.4s
Running query:   RUNNING 515.0s
Running query:   RUNNING 515.7s
Running query:   RUNNING 516.4s
Running query:   RUNNING 517.0s
Running query:   RUNNING 517.6s
Running query:   RUNNING 518.3s
Running query:   RUNNING 518.9s
Running query:   RUNNING 519.6s
Running query:   RUNNING 520.3s
Running query:   RUNNING 520.9s
Running query:   RUNNING 521.5s
Running query:   RUNNING 522.2s
Running query:   RUNNING 522.8s
Running query:   RUNNING 523.5s
Running query:   RUNNING 524.1s
Running query:   RUNNING 524.8s
Running query:   RUNNING 525.5s
Running query:   RUNNING 526.1s
Running query:   RUNNING 526.7s
Running query:   RUNNING 527.3s
Running query:   RUNNING 528.0s
Running query:   RUNNING 528.6s
Running query:   RUNNING 529.3s
Running query:   RUNNING 530.0s
Running query:   RUNNING 530.6s
Running query:   RUNNING 531.3s
Running query:   RUNNING 531.9s
Running query:   RUNNING 532.6s
Running query:   RUNNING 533.3s
Running query:   RUNNING 533.9s
Running query:   RUNNING 534.6s
Running query:   RUNNING 535.2s
Running query:   RUNNING 535.8s
Running query:   RUNNING 536.5s
Running query:   RUNNING 537.1s
Running query:   RUNNING 537.7s
Running query:   RUNNING 538.4s
Running query:   RUNNING 539.0s
Running query:   RUNNING 539.7s
Running query:   RUNNING 540.3s
Running query:   RUNNING 541.0s
Running query:   RUNNING 541.6s
Running query:   RUNNING 542.3s
Running query:   RUNNING 542.9s
Running query:   RUNNING 543.5s
Running query:   RUNNING 544.2s
Running query:   RUNNING 544.8s
Running query:   RUNNING 545.5s
Running query:   RUNNING 546.2s
Running query:   RUNNING 546.8s
Running query:   RUNNING 547.5s
Running query:   RUNNING 548.2s
Running query:   RUNNING 548.8s
Running query:   RUNNING 549.5s
Running query:   RUNNING 550.1s
Running query:   RUNNING 550.8s
Running query:   RUNNING 551.4s
Running query:   RUNNING 552.0s
Running query:   RUNNING 552.7s
Running query:   RUNNING 553.3s
Running query:   RUNNING 554.0s
Running query:   RUNNING 554.6s
Running query:   RUNNING 555.2s
Running query:   RUNNING 555.9s
Running query:   RUNNING 556.5s
Running query:   RUNNING 557.2s
Running query:   RUNNING 557.8s
Running query:   RUNNING 558.4s
Running query:   RUNNING 559.1s
Running query:   RUNNING 559.7s
Running query:   RUNNING 560.4s
Running query:   RUNNING 561.0s
Running query:   RUNNING 561.6s
Running query:   RUNNING 562.2s
Running query:   RUNNING 562.9s
Running query:   RUNNING 563.6s
Running query:   RUNNING 564.2s
Running query:   RUNNING 564.9s
Running query:   RUNNING 565.5s
Running query:   RUNNING 566.1s
Running query:   RUNNING 566.7s
Running query:   RUNNING 567.4s
Running query:   RUNNING 568.0s
Running query:   RUNNING 568.7s
Running query:   RUNNING 569.3s
Running query:   RUNNING 569.9s
Running query:   RUNNING 570.6s
Running query:   RUNNING 571.2s
Running query:   RUNNING 571.8s
Running query:   RUNNING 572.5s
Running query:   RUNNING 573.2s
Running query:   RUNNING 573.8s
Running query:   RUNNING 574.5s
Running query:   RUNNING 575.1s
Running query:   RUNNING 575.7s
Running query:   RUNNING 576.4s
Running query:   RUNNING 577.1s
Running query:   RUNNING 577.7s
Running query:   RUNNING 578.4s
Running query:   RUNNING 579.0s
Running query:   RUNNING 579.7s
Running query:   RUNNING 580.3s
Running query:   RUNNING 581.0s
Running query:   RUNNING 581.6s
Running query:   RUNNING 582.2s
Running query:   RUNNING 582.9s
Running query:   RUNNING 583.5s
Running query:   RUNNING 584.2s
Running query:   RUNNING 584.8s
Running query:   RUNNING 585.5s
Running query:   RUNNING 586.1s
Running query:   RUNNING 586.8s
Running query:   RUNNING 587.4s
Running query:   RUNNING 588.3s
Running query:   RUNNING 588.9s
Running query:   RUNNING 589.6s
Running query:   RUNNING 590.2s
Running query:   RUNNING 590.9s
Running query:   RUNNING 591.6s
Running query:   RUNNING 592.2s
Running query:   RUNNING 592.9s
Running query:   RUNNING 593.5s
Running query:   RUNNING 594.2s
Running query:   RUNNING 594.8s
Running query:   RUNNING 595.5s
Running query:   RUNNING 596.1s
Running query:   RUNNING 596.8s
Running query:   RUNNING 597.4s
Running query:   RUNNING 598.1s
Running query:   RUNNING 598.8s
Running query:   RUNNING 599.4s
Running query:   RUNNING 600.0s
Running query:   RUNNING 600.7s
Running query:   RUNNING 601.3s
Running query:   RUNNING 602.0s
Running query:   RUNNING 602.6s
Running query:   RUNNING 603.3s
Running query:   RUNNING 603.9s
Running query:   RUNNING 604.5s
Running query:   RUNNING 605.2s
Running query:   RUNNING 605.8s
Running query:   RUNNING 606.4s
Running query:   RUNNING 607.1s
Running query:   RUNNING 607.7s
Running query:   RUNNING 608.4s
Running query:   RUNNING 609.0s
Running query:   RUNNING 609.6s
Running query:   RUNNING 610.3s
Running query:   RUNNING 610.9s
Running query:   RUNNING 611.5s
Running query:   RUNNING 612.2s
Running query:   RUNNING 612.9s
Running query:   RUNNING 613.5s
Running query:   RUNNING 614.1s
Running query:   RUNNING 614.8s
Running query:   RUNNING 615.4s
Running query:   RUNNING 616.0s
Running query:   RUNNING 616.7s
Running query:   RUNNING 617.3s
Running query:   RUNNING 618.0s
Running query:   RUNNING 618.6s
Running query:   RUNNING 619.3s
Running query:   RUNNING 619.9s
Running query:   RUNNING 620.5s
Running query:   RUNNING 621.2s
Running query:   RUNNING 621.8s
Running query:   RUNNING 622.4s
Running query:   RUNNING 623.1s
Running query:   RUNNING 623.7s
Running query:   RUNNING 624.4s
Running query:   RUNNING 625.0s
Running query:   RUNNING 625.7s
Running query:   RUNNING 626.3s
Running query:   RUNNING 627.0s
Running query:   RUNNING 627.6s
Running query:   RUNNING 628.3s
Running query:   RUNNING 628.9s
Running query:   RUNNING 629.5s
Running query:   RUNNING 630.2s
Running query:   RUNNING 630.8s
Running query:   RUNNING 631.5s
Running query:   RUNNING 632.1s
Running query:   RUNNING 632.8s
Running query:   RUNNING 633.4s
Running query:   RUNNING 634.0s
Running query:   RUNNING 634.7s
Running query:   RUNNING 635.3s
Running query:   RUNNING 636.0s
Running query:   RUNNING 636.6s
Running query:   RUNNING 637.3s
Running query:   RUNNING 637.9s
Running query:   RUNNING 638.5s
Running query:   RUNNING 639.2s
Running query:   RUNNING 639.8s
Running query:   RUNNING 640.5s
Running query:   RUNNING 641.1s
Running query:   RUNNING 641.8s
Running query:   RUNNING 642.4s
Running query:   RUNNING 643.1s
Running query:   RUNNING 643.7s
Running query:   RUNNING 644.4s
Running query:   RUNNING 645.0s
Running query:   RUNNING 645.6s
Running query:   RUNNING 646.3s
Running query:   RUNNING 646.9s
Running query:   RUNNING 647.6s
Running query:   RUNNING 648.2s
Running query:   RUNNING 648.9s
Running query:   RUNNING 649.5s
Running query:   RUNNING 650.2s
Running query:   RUNNING 650.9s
Running query:   RUNNING 651.6s
```
Number of rows returned by this query: 478.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:16:30 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> calls_in_common </th> <th> identical_calls </th> <th> concordance </th>  </tr>
  <tr> <td> LP6005051-DNA_H11 </td> <td align="right"> 2197981 </td> <td align="right"> 2176158 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_A07 </td> <td align="right"> 2202307 </td> <td align="right"> 2180093 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005051-DNA_A01 </td> <td align="right"> 2201409 </td> <td align="right"> 2179276 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005144-DNA_G08 </td> <td align="right"> 2198242 </td> <td align="right"> 2176473 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005144-DNA_A10 </td> <td align="right"> 2202859 </td> <td align="right"> 2180820 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005051-DNA_G01 </td> <td align="right"> 2202919 </td> <td align="right"> 2180683 </td> <td align="right"> 0.99 </td> </tr>
   </table>


Visualizing the results:

```r
ggplot(concordanceResult) +
  geom_point(aes(x=sample_id, y=concordance)) +
  xlab("Sample") +
  ylab("Concordance") +
  ggtitle("Concordance with Genotyping Data")
```

<img src="figure/concordance-1.png" title="plot of chunk concordance" alt="plot of chunk concordance" style="display: block; margin: auto;" />

## Ethnicity Inference

TODO

## Genome Similarity

Perform a simplistic similarity check on each pair of genomes to identify any mislabled or cross-contaminated samples.

Note that this `n^2` analysis is a cluster compute job instead of a BigQuery query.




```r
require(reshape2)
require(dplyr)
ibsDataFlowFilename = '/Users/gmcinnes/data/all-genomes-ibs-2.tsv'
ReadIBSFile <- function(ibsFilename, header=FALSE, rowNames=NULL) {
  ibsData <- read.table(file=ibsFilename, header=header,
                        row.names=rowNames, stringsAsFactors=FALSE)
  return (ibsData)
}
ibsDataflowData <- ReadIBSFile(ibsDataFlowFilename)

ColumnNames <- function(ibsData) { 
  if(3 == ncol(ibsData)) {
    colnames(ibsData) <- c("sample1", "sample2", "ibsScore")
  } else {
    colnames(ibsData) <- c("sample1", "sample2", "ibsScore", "similar", "observed")
  }
}
colnames(ibsDataflowData) <- ColumnNames(ibsDataflowData)

MakeIBSDataSymmetric <- function(ibsData) {
  ibsPairsMirrored <- data.frame(sample1=ibsData$sample2,
                                 sample2=ibsData$sample1,
                                 ibsScore=ibsData$ibsScore)
  ibsData <- rbind(ibsData[,1:3], ibsPairsMirrored)
}
ibsDataflowData <- MakeIBSDataSymmetric(ibsDataflowData)

ExcludeDiagonal <- function(ibsData) {
  ibsData <- filter(ibsData, ibsData$sample1 != ibsData$sample2)
  return (ibsData)
}
ibsDataflowDataSample <- ExcludeDiagonal(ibsDataflowData)

SampleIBSMatrix <- function(ibsData, sampleSize=50) {
  individuals <- unique(ibsData$sample1)
  sample <- sample(individuals, sampleSize)
  ibsData <- subset(ibsData, ibsData$sample1 %in% sample)
  ibsData <- subset(ibsData, ibsData$sample2 %in% sample)
  return (ibsData)
}
ibsDataflowDataSubset <- SampleIBSMatrix(ibsDataflowDataSample)
```

Let's plot a subset of the data to understand the plot

```r
DrawHeatMap <- function(ibsData) {
  p <- ggplot(data=ibsData, aes(x=sample1, y=sample2)) +
    theme(axis.ticks=element_blank(), axis.text=element_blank()) +
    geom_tile(aes(fill=ibsScore), colour="white") +
    scale_fill_gradient(low="white", high="steelblue", na.value="black",
                        guide=guide_colourbar(title= "IBS Score")) +
    labs(list(title="Identity By State (IBS) Heat Map",
              x="Sample", y="Sample"))
  p
}
DrawHeatMap(ibsDataflowDataSubset)
```

<img src="figure/ibs-1-1.png" title="plot of chunk ibs-1" alt="plot of chunk ibs-1" style="display: block; margin: auto;" />

Now let's look at all the genomes


```r
DrawHeatMap(ibsDataflowDataSample)
```

<img src="figure/ibs-full-1.png" title="plot of chunk ibs-full" alt="plot of chunk ibs-full" style="display: block; margin: auto;" />

Let's take a look at the most similar genomes.
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May 15 02:16:38 2015 -->
<table border=1>
<tr> <th> sample1 </th> <th> sample2 </th> <th> ibsScore </th>  </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005692-DNA_B02 </td> <td> LP6005693-DNA_H01 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_H12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_A12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_B12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_H11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_C12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_D12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_B12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_G12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H11 </td> <td> LP6005243-DNA_F12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005693-DNA_H01 </td> <td> LP6005692-DNA_B02 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_E12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_H12 </td> <td> LP6005243-DNA_G11 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_E12 </td> <td> LP6005243-DNA_A12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_F12 </td> <td> LP6005243-DNA_D12 </td> <td align="right"> 0.10 </td> </tr>
  <tr> <td> LP6005243-DNA_G11 </td> <td> LP6005243-DNA_C12 </td> <td align="right"> 0.10 </td> </tr>
   </table>





--------------------------------------------------------
_Next_: [Part 4: Variant-Level QC](./Variant-Level-QC.md)
