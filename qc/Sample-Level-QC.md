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

*This code is used with permission by Google Genomics*
*https://github.com/googlegenomics*

# Part 3: Sample-Level QC





In Part 3 of the codelab, we perform some quality control analyses that could help to identify any problematic genomes that should be removed from the cohort before proceeding with further analysis.  The appropriate cut off thresholds will depend upon the input dataset and/or other factors.

By default this codelab runs upon the Illumina Platinum Genomes Variants. Update the table and change the source of sample information here if you wish to run the queries against a different dataset.

```r
tableReplacement <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs",
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
      [va_aaa_pilot_data.all_genomes_gvcfs])
    GROUP BY
      call.call_set_name)) AS g
  CROSS JOIN (
    SELECT 
      COUNT(Chr) AS count
    FROM 
      [google.com:biggene:test.hg19]) AS hg19
ORDER BY
  missingness DESC
```
Number of rows returned by this query: 480.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Tue May 12 23:21:14 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> missingness </th>  </tr>
  <tr> <td> LP6005243-DNA_A08 </td> <td align="right"> 0.55 </td> </tr>
  <tr> <td> LP6005038-DNA_C02 </td> <td align="right"> 0.06 </td> </tr>
  <tr> <td> LP6005692-DNA_D09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005243-DNA_C01 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005243-DNA_D01 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A01 </td> <td align="right"> 0.05 </td> </tr>
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
            [va_aaa_pilot_data.all_genomes_gvcfs]
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
Number of rows returned by this query: 480.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Tue May 12 23:21:17 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> private_variant_count </th>  </tr>
  <tr> <td> LP6005051-DNA_D09 </td> <td align="right"> 498164 </td> </tr>
  <tr> <td> LP6005692-DNA_D05 </td> <td align="right"> 403557 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td align="right"> 170457 </td> </tr>
  <tr> <td> LP6005692-DNA_H01 </td> <td align="right"> 137485 </td> </tr>
  <tr> <td> LP6005144-DNA_D03 </td> <td align="right"> 136157 </td> </tr>
  <tr> <td> LP6005243-DNA_G07 </td> <td align="right"> 81819 </td> </tr>
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
```
Number of rows returned by this query: 479.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Tue May 12 23:21:20 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 25389354 </td> <td align="right"> 29665850.23 </td> <td align="right"> 27556858 </td> <td align="right"> 2.03 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 25228070 </td> <td align="right"> 29315956.80 </td> <td align="right"> 27395111 </td> <td align="right"> 2.13 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 25178694 </td> <td align="right"> 29370309.77 </td> <td align="right"> 27429303 </td> <td align="right"> 2.16 </td> </tr>
  <tr> <td> LP6005038-DNA_A04 </td> <td align="right"> 25291057 </td> <td align="right"> 29449908.26 </td> <td align="right"> 27445069 </td> <td align="right"> 2.07 </td> </tr>
  <tr> <td> LP6005038-DNA_A05 </td> <td align="right"> 25102231 </td> <td align="right"> 29255773.98 </td> <td align="right"> 27347464 </td> <td align="right"> 2.18 </td> </tr>
  <tr> <td> LP6005038-DNA_A06 </td> <td align="right"> 25227077 </td> <td align="right"> 29409374.64 </td> <td align="right"> 27430575 </td> <td align="right"> 2.11 </td> </tr>
   </table>

And visualizing the results:

```r
limits <- c(min(result$O_HOM, result$E_HOM),
            max(result$O_HOM, result$E_HOM))

ggplot(result) +
  geom_point(aes(x=O_HOM, y=E_HOM, label=call_call_set_name), alpha=1/1.5) +
  geom_abline(color="darkslateblue") +
  scale_x_continuous(limits=limits) + 
  scale_y_continuous(limits=limits) +
  xlab("Observed Homozygous Variants") +
  ylab("Expected Homozygous Variants") +
  ggtitle("Homozygosity")
```

<img src="figure/homozygosity-1.png" title="plot of chunk homozygosity" alt="plot of chunk homozygosity" style="display: block; margin: auto;" />

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

<img src="figure/homozygosity-labeled-1.png" title="plot of chunk homozygosity-labeled" alt="plot of chunk homozygosity-labeled" style="display: block; margin: auto;" />


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
<!-- Tue May 12 23:21:23 2015 -->
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
          [va_aaa_pilot_data.all_genomes_gvcfs]
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
Running query:   RUNNING  2.7s
Running query:   RUNNING  3.3s
Running query:   RUNNING  4.0s
Running query:   RUNNING  4.6s
Running query:   RUNNING  5.3s
Running query:   RUNNING  5.9s
Running query:   RUNNING  6.5s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.8s
Running query:   RUNNING  8.4s
Running query:   RUNNING  9.0s
Running query:   RUNNING  9.6s
Running query:   RUNNING 10.3s
Running query:   RUNNING 10.9s
Running query:   RUNNING 11.5s
Running query:   RUNNING 12.1s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.4s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.6s
Running query:   RUNNING 15.2s
Running query:   RUNNING 15.9s
Running query:   RUNNING 16.5s
Running query:   RUNNING 17.1s
Running query:   RUNNING 17.7s
Running query:   RUNNING 18.4s
Running query:   RUNNING 19.0s
Running query:   RUNNING 19.7s
Running query:   RUNNING 20.3s
Running query:   RUNNING 20.9s
Running query:   RUNNING 21.5s
Running query:   RUNNING 22.1s
Running query:   RUNNING 22.8s
Running query:   RUNNING 23.4s
Running query:   RUNNING 24.0s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.2s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.5s
Running query:   RUNNING 27.1s
Running query:   RUNNING 27.8s
Running query:   RUNNING 28.4s
Running query:   RUNNING 29.0s
Running query:   RUNNING 29.6s
Running query:   RUNNING 30.2s
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
Running query:   RUNNING 37.1s
Running query:   RUNNING 37.7s
Running query:   RUNNING 38.4s
Running query:   RUNNING 39.0s
Running query:   RUNNING 39.6s
Running query:   RUNNING 40.2s
Running query:   RUNNING 40.9s
Running query:   RUNNING 41.5s
Running query:   RUNNING 42.1s
Running query:   RUNNING 42.8s
Running query:   RUNNING 43.4s
Running query:   RUNNING 44.0s
Running query:   RUNNING 44.6s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.9s
Running query:   RUNNING 46.5s
Running query:   RUNNING 47.1s
Running query:   RUNNING 47.7s
Running query:   RUNNING 48.3s
Running query:   RUNNING 48.9s
Running query:   RUNNING 49.6s
Running query:   RUNNING 50.2s
Running query:   RUNNING 50.8s
Running query:   RUNNING 51.5s
Running query:   RUNNING 52.1s
Running query:   RUNNING 52.7s
Running query:   RUNNING 53.3s
Running query:   RUNNING 54.0s
Running query:   RUNNING 54.6s
Running query:   RUNNING 55.2s
Running query:   RUNNING 55.8s
Running query:   RUNNING 56.5s
Running query:   RUNNING 57.1s
Running query:   RUNNING 57.7s
Running query:   RUNNING 58.3s
Running query:   RUNNING 58.9s
Running query:   RUNNING 59.6s
Running query:   RUNNING 60.2s
Running query:   RUNNING 60.8s
Running query:   RUNNING 61.4s
Running query:   RUNNING 62.0s
Running query:   RUNNING 62.7s
Running query:   RUNNING 63.3s
Running query:   RUNNING 64.0s
Running query:   RUNNING 64.6s
Running query:   RUNNING 65.2s
Running query:   RUNNING 65.9s
Running query:   RUNNING 66.5s
Running query:   RUNNING 67.2s
Running query:   RUNNING 67.8s
Running query:   RUNNING 68.4s
Running query:   RUNNING 69.0s
Running query:   RUNNING 69.6s
Running query:   RUNNING 70.3s
Running query:   RUNNING 70.9s
Running query:   RUNNING 71.5s
Running query:   RUNNING 72.1s
Running query:   RUNNING 72.7s
Running query:   RUNNING 73.4s
Running query:   RUNNING 74.1s
Running query:   RUNNING 74.8s
Running query:   RUNNING 75.4s
Running query:   RUNNING 76.0s
Running query:   RUNNING 76.7s
Running query:   RUNNING 77.3s
Running query:   RUNNING 77.9s
Running query:   RUNNING 78.5s
Running query:   RUNNING 79.2s
Running query:   RUNNING 79.8s
Running query:   RUNNING 80.4s
Running query:   RUNNING 81.0s
Running query:   RUNNING 81.7s
Running query:   RUNNING 82.3s
Running query:   RUNNING 82.9s
Running query:   RUNNING 83.5s
Running query:   RUNNING 84.1s
Running query:   RUNNING 84.7s
Running query:   RUNNING 85.4s
Running query:   RUNNING 86.0s
Running query:   RUNNING 86.7s
Running query:   RUNNING 87.3s
Running query:   RUNNING 87.9s
Running query:   RUNNING 88.6s
Running query:   RUNNING 89.2s
Running query:   RUNNING 89.9s
Running query:   RUNNING 90.5s
Running query:   RUNNING 91.2s
Running query:   RUNNING 91.8s
Running query:   RUNNING 92.5s
Running query:   RUNNING 93.1s
Running query:   RUNNING 93.7s
Running query:   RUNNING 94.3s
Running query:   RUNNING 94.9s
Running query:   RUNNING 95.6s
Running query:   RUNNING 96.2s
Running query:   RUNNING 96.8s
Running query:   RUNNING 97.5s
Running query:   RUNNING 98.1s
Running query:   RUNNING 98.7s
Running query:   RUNNING 99.4s
Running query:   RUNNING 100.0s
Running query:   RUNNING 100.6s
Running query:   RUNNING 101.2s
Running query:   RUNNING 101.9s
Running query:   RUNNING 102.6s
Running query:   RUNNING 103.2s
Running query:   RUNNING 103.8s
Running query:   RUNNING 104.4s
Running query:   RUNNING 105.0s
Running query:   RUNNING 105.7s
Running query:   RUNNING 106.3s
Running query:   RUNNING 106.9s
Running query:   RUNNING 107.6s
Running query:   RUNNING 108.2s
Running query:   RUNNING 108.8s
Running query:   RUNNING 109.5s
Running query:   RUNNING 110.1s
Running query:   RUNNING 110.9s
Running query:   RUNNING 111.5s
Running query:   RUNNING 112.1s
Running query:   RUNNING 112.8s
Running query:   RUNNING 113.4s
Running query:   RUNNING 114.0s
Running query:   RUNNING 114.6s
Running query:   RUNNING 115.3s
Running query:   RUNNING 116.4s
Running query:   RUNNING 117.1s
Running query:   RUNNING 117.7s
Running query:   RUNNING 118.3s
Running query:   RUNNING 118.9s
Running query:   RUNNING 119.6s
Running query:   RUNNING 120.2s
Running query:   RUNNING 120.8s
Running query:   RUNNING 121.4s
Running query:   RUNNING 122.1s
Running query:   RUNNING 122.7s
Running query:   RUNNING 123.3s
Running query:   RUNNING 123.9s
Running query:   RUNNING 124.5s
Running query:   RUNNING 125.1s
Running query:   RUNNING 125.8s
Running query:   RUNNING 126.4s
Running query:   RUNNING 127.0s
Running query:   RUNNING 127.6s
Running query:   RUNNING 128.3s
Running query:   RUNNING 128.9s
Running query:   RUNNING 129.5s
Running query:   RUNNING 130.2s
Running query:   RUNNING 130.8s
Running query:   RUNNING 131.4s
Running query:   RUNNING 132.0s
Running query:   RUNNING 133.3s
Running query:   RUNNING 134.1s
Running query:   RUNNING 134.7s
Running query:   RUNNING 135.3s
Running query:   RUNNING 136.0s
Running query:   RUNNING 136.6s
Running query:   RUNNING 137.2s
Running query:   RUNNING 137.9s
Running query:   RUNNING 138.5s
Running query:   RUNNING 139.1s
Running query:   RUNNING 139.7s
Running query:   RUNNING 140.4s
Running query:   RUNNING 141.0s
Running query:   RUNNING 141.6s
Running query:   RUNNING 142.2s
Running query:   RUNNING 142.9s
Running query:   RUNNING 143.5s
Running query:   RUNNING 144.1s
Running query:   RUNNING 144.8s
Running query:   RUNNING 145.4s
Running query:   RUNNING 146.0s
Running query:   RUNNING 146.6s
Running query:   RUNNING 147.2s
Running query:   RUNNING 147.9s
Running query:   RUNNING 148.5s
Running query:   RUNNING 149.1s
Running query:   RUNNING 149.8s
Running query:   RUNNING 150.4s
Running query:   RUNNING 151.1s
Running query:   RUNNING 151.7s
Running query:   RUNNING 152.3s
Running query:   RUNNING 152.9s
Running query:   RUNNING 153.5s
Running query:   RUNNING 154.1s
Running query:   RUNNING 154.8s
Running query:   RUNNING 155.4s
Running query:   RUNNING 156.0s
Running query:   RUNNING 156.6s
Running query:   RUNNING 157.3s
Running query:   RUNNING 157.9s
Running query:   RUNNING 158.5s
Running query:   RUNNING 159.1s
Running query:   RUNNING 159.8s
Running query:   RUNNING 160.4s
Running query:   RUNNING 161.0s
Running query:   RUNNING 161.6s
Running query:   RUNNING 162.2s
Running query:   RUNNING 162.9s
Running query:   RUNNING 163.5s
Running query:   RUNNING 164.1s
Running query:   RUNNING 164.8s
Running query:   RUNNING 165.4s
Running query:   RUNNING 166.0s
Running query:   RUNNING 166.6s
Running query:   RUNNING 167.2s
Running query:   RUNNING 167.9s
Running query:   RUNNING 168.5s
Running query:   RUNNING 169.1s
Running query:   RUNNING 169.7s
Running query:   RUNNING 170.4s
Running query:   RUNNING 171.0s
Running query:   RUNNING 171.6s
Running query:   RUNNING 172.3s
Running query:   RUNNING 172.9s
Running query:   RUNNING 173.5s
Running query:   RUNNING 174.1s
Running query:   RUNNING 174.8s
Running query:   RUNNING 175.4s
Running query:   RUNNING 176.0s
Running query:   RUNNING 176.6s
Running query:   RUNNING 177.2s
Running query:   RUNNING 177.9s
Running query:   RUNNING 178.5s
Running query:   RUNNING 179.2s
Running query:   RUNNING 179.8s
Running query:   RUNNING 180.4s
Running query:   RUNNING 181.0s
Running query:   RUNNING 181.7s
Running query:   RUNNING 182.3s
Running query:   RUNNING 182.9s
Running query:   RUNNING 183.5s
Running query:   RUNNING 184.2s
Running query:   RUNNING 184.8s
Running query:   RUNNING 185.4s
Running query:   RUNNING 186.0s
Running query:   RUNNING 186.7s
Running query:   RUNNING 187.3s
Running query:   RUNNING 187.9s
Running query:   RUNNING 188.5s
Running query:   RUNNING 189.1s
Running query:   RUNNING 189.8s
Running query:   RUNNING 190.4s
Running query:   RUNNING 191.1s
Running query:   RUNNING 191.7s
Running query:   RUNNING 192.4s
Running query:   RUNNING 193.0s
Running query:   RUNNING 193.6s
Running query:   RUNNING 194.3s
Running query:   RUNNING 194.9s
Running query:   RUNNING 195.5s
Running query:   RUNNING 196.2s
Running query:   RUNNING 196.8s
Running query:   RUNNING 197.5s
Running query:   RUNNING 198.1s
Running query:   RUNNING 198.7s
Running query:   RUNNING 199.3s
Running query:   RUNNING 199.9s
Running query:   RUNNING 200.6s
Running query:   RUNNING 201.2s
Running query:   RUNNING 201.9s
Running query:   RUNNING 202.5s
Running query:   RUNNING 203.1s
Running query:   RUNNING 203.8s
Running query:   RUNNING 204.4s
Running query:   RUNNING 205.0s
Running query:   RUNNING 205.6s
Running query:   RUNNING 206.3s
Running query:   RUNNING 206.9s
Running query:   RUNNING 207.8s
Running query:   RUNNING 208.4s
Running query:   RUNNING 209.0s
Running query:   RUNNING 209.6s
Running query:   RUNNING 210.3s
Running query:   RUNNING 210.9s
Running query:   RUNNING 211.5s
Running query:   RUNNING 212.1s
Running query:   RUNNING 212.8s
Running query:   RUNNING 213.4s
Running query:   RUNNING 214.0s
Running query:   RUNNING 214.7s
Running query:   RUNNING 215.3s
Running query:   RUNNING 215.9s
Running query:   RUNNING 216.6s
Running query:   RUNNING 217.2s
Running query:   RUNNING 217.8s
Running query:   RUNNING 218.4s
Running query:   RUNNING 219.1s
Running query:   RUNNING 219.8s
Running query:   RUNNING 220.4s
Running query:   RUNNING 221.0s
Running query:   RUNNING 221.6s
Running query:   RUNNING 222.3s
Running query:   RUNNING 222.9s
Running query:   RUNNING 223.6s
Running query:   RUNNING 224.3s
Running query:   RUNNING 224.9s
Running query:   RUNNING 225.6s
Running query:   RUNNING 226.3s
Running query:   RUNNING 226.9s
Running query:   RUNNING 227.6s
Running query:   RUNNING 228.2s
Running query:   RUNNING 228.9s
Running query:   RUNNING 229.5s
Running query:   RUNNING 230.2s
Running query:   RUNNING 230.8s
Running query:   RUNNING 231.5s
Running query:   RUNNING 232.1s
Running query:   RUNNING 232.7s
Running query:   RUNNING 233.4s
Running query:   RUNNING 234.0s
Running query:   RUNNING 234.6s
Running query:   RUNNING 235.3s
Running query:   RUNNING 235.9s
Running query:   RUNNING 236.6s
Running query:   RUNNING 237.2s
Running query:   RUNNING 237.8s
Running query:   RUNNING 238.4s
Running query:   RUNNING 239.1s
Running query:   RUNNING 239.7s
Running query:   RUNNING 240.4s
Running query:   RUNNING 241.0s
Running query:   RUNNING 241.7s
Running query:   RUNNING 242.3s
Running query:   RUNNING 243.0s
Running query:   RUNNING 243.6s
Running query:   RUNNING 244.3s
Running query:   RUNNING 244.9s
Running query:   RUNNING 245.7s
Running query:   RUNNING 246.3s
Running query:   RUNNING 247.0s
Running query:   RUNNING 247.7s
Running query:   RUNNING 248.3s
Running query:   RUNNING 249.0s
Running query:   RUNNING 249.6s
Running query:   RUNNING 250.3s
Running query:   RUNNING 250.9s
Running query:   RUNNING 251.6s
Running query:   RUNNING 252.3s
Running query:   RUNNING 252.9s
Running query:   RUNNING 253.6s
Running query:   RUNNING 254.2s
Running query:   RUNNING 254.9s
Running query:   RUNNING 255.5s
Running query:   RUNNING 256.2s
Running query:   RUNNING 256.8s
Running query:   RUNNING 257.5s
Running query:   RUNNING 258.1s
Running query:   RUNNING 258.8s
Running query:   RUNNING 259.5s
Running query:   RUNNING 260.1s
Running query:   RUNNING 260.7s
Running query:   RUNNING 261.4s
Running query:   RUNNING 262.1s
Running query:   RUNNING 262.7s
Running query:   RUNNING 263.4s
Running query:   RUNNING 264.0s
Running query:   RUNNING 264.7s
Running query:   RUNNING 265.3s
Running query:   RUNNING 266.0s
Running query:   RUNNING 266.7s
Running query:   RUNNING 267.3s
Running query:   RUNNING 267.9s
Running query:   RUNNING 268.6s
Running query:   RUNNING 269.2s
Running query:   RUNNING 269.9s
Running query:   RUNNING 270.5s
Running query:   RUNNING 271.2s
Running query:   RUNNING 271.8s
Running query:   RUNNING 272.5s
Running query:   RUNNING 273.1s
Running query:   RUNNING 273.8s
Running query:   RUNNING 274.4s
Running query:   RUNNING 275.0s
Running query:   RUNNING 275.6s
Running query:   RUNNING 276.3s
Running query:   RUNNING 276.9s
Running query:   RUNNING 277.5s
Running query:   RUNNING 278.2s
Running query:   RUNNING 278.8s
Running query:   RUNNING 279.4s
Running query:   RUNNING 280.1s
Running query:   RUNNING 280.7s
Running query:   RUNNING 281.3s
Running query:   RUNNING 281.9s
Running query:   RUNNING 282.6s
Running query:   RUNNING 283.3s
Running query:   RUNNING 283.9s
Running query:   RUNNING 284.5s
Running query:   RUNNING 285.2s
Running query:   RUNNING 285.8s
Running query:   RUNNING 286.4s
Running query:   RUNNING 287.1s
Running query:   RUNNING 287.7s
Running query:   RUNNING 288.4s
Running query:   RUNNING 289.0s
Running query:   RUNNING 289.6s
Running query:   RUNNING 290.3s
Running query:   RUNNING 290.9s
Running query:   RUNNING 291.5s
Running query:   RUNNING 292.2s
Running query:   RUNNING 292.8s
Running query:   RUNNING 293.5s
Running query:   RUNNING 294.1s
Running query:   RUNNING 294.7s
Running query:   RUNNING 295.3s
Running query:   RUNNING 296.0s
Running query:   RUNNING 296.6s
Running query:   RUNNING 297.3s
Running query:   RUNNING 297.9s
Running query:   RUNNING 298.5s
Running query:   RUNNING 299.2s
Running query:   RUNNING 299.8s
Running query:   RUNNING 300.4s
Running query:   RUNNING 301.1s
Running query:   RUNNING 301.8s
Running query:   RUNNING 302.4s
Running query:   RUNNING 303.0s
Running query:   RUNNING 303.7s
Running query:   RUNNING 304.3s
Running query:   RUNNING 304.9s
Running query:   RUNNING 305.6s
Running query:   RUNNING 306.2s
Running query:   RUNNING 306.8s
Running query:   RUNNING 307.4s
Running query:   RUNNING 308.1s
Running query:   RUNNING 308.7s
Running query:   RUNNING 309.3s
Running query:   RUNNING 310.0s
Running query:   RUNNING 310.6s
Running query:   RUNNING 311.3s
Running query:   RUNNING 311.9s
Running query:   RUNNING 312.5s
Running query:   RUNNING 313.2s
Running query:   RUNNING 313.8s
Running query:   RUNNING 314.4s
Running query:   RUNNING 315.1s
Running query:   RUNNING 315.7s
Running query:   RUNNING 316.3s
Running query:   RUNNING 317.0s
Running query:   RUNNING 317.6s
Running query:   RUNNING 318.2s
Running query:   RUNNING 318.9s
Running query:   RUNNING 319.5s
Running query:   RUNNING 320.2s
Running query:   RUNNING 320.8s
Running query:   RUNNING 321.4s
Running query:   RUNNING 322.1s
Running query:   RUNNING 322.7s
Running query:   RUNNING 323.3s
Running query:   RUNNING 324.0s
Running query:   RUNNING 324.6s
Running query:   RUNNING 325.2s
Running query:   RUNNING 325.9s
Running query:   RUNNING 326.6s
Running query:   RUNNING 327.2s
Running query:   RUNNING 327.8s
Running query:   RUNNING 328.4s
Running query:   RUNNING 329.1s
Running query:   RUNNING 329.7s
Running query:   RUNNING 330.3s
Running query:   RUNNING 330.9s
Running query:   RUNNING 331.6s
Running query:   RUNNING 332.2s
Running query:   RUNNING 332.8s
Running query:   RUNNING 333.5s
Running query:   RUNNING 334.1s
Running query:   RUNNING 334.8s
Running query:   RUNNING 335.4s
Running query:   RUNNING 336.1s
Running query:   RUNNING 336.7s
Running query:   RUNNING 337.3s
Running query:   RUNNING 338.0s
Running query:   RUNNING 338.6s
Running query:   RUNNING 339.3s
Running query:   RUNNING 339.9s
Running query:   RUNNING 340.6s
```
Number of rows returned by this query: 478.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Tue May 12 23:27:08 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> calls_in_common </th> <th> identical_calls </th> <th> concordance </th>  </tr>
  <tr> <td> LP6005038-DNA_H05 </td> <td align="right"> 2200555 </td> <td align="right"> 2178550 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005692-DNA_C01 </td> <td align="right"> 2177148 </td> <td align="right"> 2140591 </td> <td align="right"> 0.98 </td> </tr>
  <tr> <td> LP6005051-DNA_D02 </td> <td align="right"> 2200793 </td> <td align="right"> 2178731 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_G12 </td> <td align="right"> 2201643 </td> <td align="right"> 2179570 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005144-DNA_E09 </td> <td align="right"> 2200995 </td> <td align="right"> 2178761 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005051-DNA_A02 </td> <td align="right"> 2199365 </td> <td align="right"> 2177646 </td> <td align="right"> 0.99 </td> </tr>
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
ibsDataflowDataSample <- SampleIBSMatrix(ibsDataflowDataSample)
```


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
DrawHeatMap(ibsDataflowDataSample)
```

<img src="figure/ibs-1-1.png" title="plot of chunk ibs-1" alt="plot of chunk ibs-1" style="display: block; margin: auto;" />


Let's take a look at the most similar genomes.

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Tue May 12 23:27:11 2015 -->
<table border=1>
<tr> <th> sample1 </th> <th> sample2 </th> <th> ibsScore </th>  </tr>
  <tr> <td> LP6005038-DNA_F08 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_A04 </td> <td> LP6005692-DNA_E08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_D04 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005692-DNA_E08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_E08 </td> <td> LP6005692-DNA_H12 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_E08 </td> <td> LP6005693-DNA_D01 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_G07 </td> <td> LP6005692-DNA_H12 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_A04 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005144-DNA_A09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_G09 </td> <td> LP6005692-DNA_H11 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_G07 </td> <td> LP6005692-DNA_E08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005693-DNA_D01 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H11 </td> <td> LP6005693-DNA_D01 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_C09 </td> <td> LP6005692-DNA_H12 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005692-DNA_H12 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005243-DNA_A05 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005051-DNA_A02 </td> <td> LP6005692-DNA_H11 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_B05 </td> <td> LP6005692-DNA_H11 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_C09 </td> <td> LP6005038-DNA_G07 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_G07 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005692-DNA_H11 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_A01 </td> <td> LP6005692-DNA_H12 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005038-DNA_F08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_E08 </td> <td> LP6005051-DNA_A04 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005038-DNA_D04 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_E08 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H12 </td> <td> LP6005692-DNA_E08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005693-DNA_D01 </td> <td> LP6005692-DNA_E08 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H12 </td> <td> LP6005038-DNA_G07 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005051-DNA_A04 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A09 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H11 </td> <td> LP6005051-DNA_G09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_E08 </td> <td> LP6005038-DNA_G07 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005693-DNA_D01 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005693-DNA_D01 </td> <td> LP6005692-DNA_H11 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H12 </td> <td> LP6005038-DNA_C09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H12 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005243-DNA_A05 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H11 </td> <td> LP6005051-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H11 </td> <td> LP6005692-DNA_B05 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005038-DNA_G07 </td> <td> LP6005038-DNA_C09 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005144-DNA_A02 </td> <td> LP6005038-DNA_G07 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H11 </td> <td> LP6005144-DNA_A02 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> LP6005692-DNA_H12 </td> <td> LP6005692-DNA_A01 </td> <td align="right"> 0.05 </td> </tr>
   </table>





--------------------------------------------------------
_Next_: [Part 4: Variant-Level QC](./Variant-Level-QC.md)
