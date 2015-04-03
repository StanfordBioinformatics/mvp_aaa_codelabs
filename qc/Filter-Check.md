# Filtering Validation

Here we show that the tables created with the calls filtered by [gvcf-mapper.py](https://github.com/StanfordBioinformatics/googva/blob/master/gvcf-mapper.py) contain all the calls that meet our filtering criteria and nothing extra.  For details on how this table was prepared see [here](./Data-Preparation-GC.md).

* [Setup](#setup)
* [Genomes in dataset](#genomes-in-dataset)
* [Variant filters](#variant-filters)
* [Reference filteres](#reference-filters)
* [Variant Details](#variant-details)

## Setup






```r
# By default this codelab runs upon the Illumina Platinum Genomes Variants.  
# Change the table here if you wish to run these queries against your own data.
queryReplacements <- list("_THE_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2")

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpCGIOnlyDataset.R")
```


## Genomes in dataset
Let's first make sure all of our samples made it into BigQuery.


```r
query <- "./sql/sample-names.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
SELECT
  call.call_set_name AS sample_name,
FROM 
  [va_aaa_pilot_data.5_genome_test_gvcfs_2]
GROUP BY
  sample_name,
ORDER BY
  sample_name,
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:19 2015 -->
<table border=1>
<tr> <th> sample_name </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> </tr>
   </table>

Looks good.  All of our samples are there.
  
## Variant filters

```r
query <- "./sql/variant-filters.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
SELECT
  call.call_set_name AS sample_name,
  call.FILTER AS filter,
  COUNT(call.FILTER) AS count
FROM (
  SELECT
    call.call_set_name,
    call.FILTER,
    alternate_bases,
    GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
  FROM (
    FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)))
WHERE 
  genotype != '-1,-1'
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  sample_name,
  filter
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:21 2015 -->
<table border=1>
<tr> <th> sample_name </th> <th> filter </th> <th> count </th>  </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td> PASS </td> <td align="right"> 4147334 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td> PASS </td> <td align="right"> 4261151 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> PASS </td> <td align="right"> 4147994 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td> PASS </td> <td align="right"> 4199099 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td> PASS </td> <td align="right"> 4199487 </td> </tr>
   </table>

## Ref filters
Let's look at the reference calls and make sure they all meet the filtering requirements.

* QUAL > 30
* MQ > 30
* MQ0 < 4


```r
query <- "./sql/ref-quality.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
SELECT
  call.call_set_name AS sample_name,
  MIN(call.QUAL) AS min_QUAL,
  MIN(MQ) AS min_MQ,
  MAX(MQ0) AS max_MQ0,
  COUNT(call.call_set_name) AS count
FROM (
  SELECT
    call.call_set_name,
    alternate_bases,
    call.QUAL,
    MQ,
    MQ0,
    GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
  FROM (
    FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)))
WHERE 
  genotype != '-1,-1'
OMIT RECORD IF EVERY(alternate_bases IS NOT NULL)
GROUP BY
  sample_name,
ORDER BY
  sample_name,
  count DESC
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:23 2015 -->
<table border=1>
<tr> <th> sample_name </th> <th> min_QUAL </th> <th> min_MQ </th> <th> max_MQ0 </th> <th> count </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 30.24 </td> <td align="right"> 30.00 </td> <td align="right">   3 </td> <td align="right"> 5087052 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 30.24 </td> <td align="right"> 10.12 </td> <td align="right">  19 </td> <td align="right"> 5152377 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 30.23 </td> <td align="right"> 9.55 </td> <td align="right">  30 </td> <td align="right"> 5030738 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 30.24 </td> <td align="right"> 7.60 </td> <td align="right">  30 </td> <td align="right"> 5048016 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 30.24 </td> <td align="right"> 6.24 </td> <td align="right">  33 </td> <td align="right"> 5063135 </td> </tr>
   </table>

We can see from the table that each of our samples has a QUAL score greater than 30, but the MQ an MQ0 columns have values outside our desired range.  This is because when variants are imported into variant store in Google the value fields are merged for identical sites.  Since we have no-calls (calls outside our desired quality range) in our dataset good calls and bad calls are being merged together.  This is why it is important to filter our data prior to importing to variant store. We can see that the first sample passes for each metric.  Subsequent samples have calls that are merged with no-calls from the first sample leading to values ourside of our desired range.  All samples pass for QUAL because that metric is stored for each call within 'call'.  We can be sure that the filtering is working properly based on the first sample.

## Variant details
all of these are messed up.  show brca1 first (which looks good) then show whole genome.  why are they different?????
Check that count of 1,2 genotypes does not match difference, pretty sure it doesn't but check

### SNPs
#### BRCA1
##### Count


```r
query <- "./sql/snp-count.sql"
where <- list("#_AND_"="AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements,where))
```

```
# overly complicated but does the job

SELECT
  COUNT(mutation) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
      OMIT RECORD if every(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
    )
  WHERE
    genotype != '-1,-1'
    AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:24 2015 -->
<table border=1>
<tr> <th> count </th>  </tr>
  <tr> <td align="right"> 138 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/mutation_counts_brca1.csv'
expectedResult = read.csv(file)
bcftoolsCount = expectedResult[expectedResult$type=='SNPs',]
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:24 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> SNPs </td> <td align="right"> 138 </td> </tr>
   </table>

##### Substitution types


```r
query <- "./sql/unique-snp-positions.sql"
where <- list("#_AND_"="AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, where))
```

```
# Count the number of positions for each type of mutation

SELECT
  mutation,
  COUNT(mutation) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD if every(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T'))
  WHERE 
    genotype != '-1,-1'
    AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype)
GROUP EACH BY mutation
ORDER BY mutation
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:27 2015 -->
<table border=1>
<tr> <th> mutation </th> <th> count </th>  </tr>
  <tr> <td> A-&gt;C </td> <td align="right">   8 </td> </tr>
  <tr> <td> A-&gt;G </td> <td align="right">  22 </td> </tr>
  <tr> <td> A-&gt;T </td> <td align="right">   6 </td> </tr>
  <tr> <td> C-&gt;A </td> <td align="right">   3 </td> </tr>
  <tr> <td> C-&gt;G </td> <td align="right">   1 </td> </tr>
  <tr> <td> C-&gt;T </td> <td align="right">  23 </td> </tr>
  <tr> <td> G-&gt;A </td> <td align="right">  22 </td> </tr>
  <tr> <td> G-&gt;C </td> <td align="right">   7 </td> </tr>
  <tr> <td> G-&gt;T </td> <td align="right">   3 </td> </tr>
  <tr> <td> T-&gt;A </td> <td align="right">   4 </td> </tr>
  <tr> <td> T-&gt;C </td> <td align="right">  33 </td> </tr>
  <tr> <td> T-&gt;G </td> <td align="right">   6 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/snp_types_brca1.csv'
expectedResult = read.csv(file)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:27 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> A&gt;C </td> <td align="right">   8 </td> </tr>
  <tr> <td> A&gt;G </td> <td align="right">  22 </td> </tr>
  <tr> <td> A&gt;T </td> <td align="right">   6 </td> </tr>
  <tr> <td> C&gt;A </td> <td align="right">   3 </td> </tr>
  <tr> <td> C&gt;G </td> <td align="right">   1 </td> </tr>
  <tr> <td> C&gt;T </td> <td align="right">  23 </td> </tr>
  <tr> <td> G&gt;A </td> <td align="right">  22 </td> </tr>
  <tr> <td> G&gt;C </td> <td align="right">   7 </td> </tr>
  <tr> <td> G&gt;T </td> <td align="right">   3 </td> </tr>
  <tr> <td> T&gt;A </td> <td align="right">   4 </td> </tr>
  <tr> <td> T&gt;C </td> <td align="right">  33 </td> </tr>
  <tr> <td> T&gt;G </td> <td align="right">   6 </td> </tr>
   </table>



#### Full Genomes
##### Count


```r
query <- "./sql/snp-count.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
# overly complicated but does the job

SELECT
  COUNT(mutation) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
      OMIT RECORD if every(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
    )
  WHERE
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:29 2015 -->
<table border=1>
<tr> <th> count </th>  </tr>
  <tr> <td align="right"> 8870189 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/mutation_counts.csv'
expectedResult = read.csv(file)
bcftoolsCount = expectedResult[expectedResult$type=='SNPs',]
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:29 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> SNPs </td> <td align="right"> 6623215 </td> </tr>
   </table>

##### Substitution types


```r
query <- "./sql/unique-snp-positions.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of positions for each type of mutation

SELECT
  mutation,
  COUNT(mutation) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD if every(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T'))
  WHERE 
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype)
GROUP EACH BY mutation
ORDER BY mutation
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:31 2015 -->
<table border=1>
<tr> <th> mutation </th> <th> count </th>  </tr>
  <tr> <td> A-&gt;C </td> <td align="right"> 355520 </td> </tr>
  <tr> <td> A-&gt;G </td> <td align="right"> 1414192 </td> </tr>
  <tr> <td> A-&gt;T </td> <td align="right"> 311282 </td> </tr>
  <tr> <td> C-&gt;A </td> <td align="right"> 382908 </td> </tr>
  <tr> <td> C-&gt;G </td> <td align="right"> 384540 </td> </tr>
  <tr> <td> C-&gt;T </td> <td align="right"> 1583313 </td> </tr>
  <tr> <td> G-&gt;A </td> <td align="right"> 1584918 </td> </tr>
  <tr> <td> G-&gt;C </td> <td align="right"> 385076 </td> </tr>
  <tr> <td> G-&gt;T </td> <td align="right"> 384424 </td> </tr>
  <tr> <td> T-&gt;A </td> <td align="right"> 311115 </td> </tr>
  <tr> <td> T-&gt;C </td> <td align="right"> 1417409 </td> </tr>
  <tr> <td> T-&gt;G </td> <td align="right"> 355492 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/snp_types.csv'
expectedResult = read.csv(file)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:31 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> A&gt;C </td> <td align="right"> 262450 </td> </tr>
  <tr> <td> A&gt;G </td> <td align="right"> 1039552 </td> </tr>
  <tr> <td> A&gt;T </td> <td align="right"> 232973 </td> </tr>
  <tr> <td> C&gt;A </td> <td align="right"> 288685 </td> </tr>
  <tr> <td> C&gt;G </td> <td align="right"> 287387 </td> </tr>
  <tr> <td> C&gt;T </td> <td align="right"> 1201255 </td> </tr>
  <tr> <td> G&gt;A </td> <td align="right"> 1203262 </td> </tr>
  <tr> <td> G&gt;C </td> <td align="right"> 287456 </td> </tr>
  <tr> <td> G&gt;T </td> <td align="right"> 289835 </td> </tr>
  <tr> <td> T&gt;A </td> <td align="right"> 232721 </td> </tr>
  <tr> <td> T&gt;C </td> <td align="right"> 1041834 </td> </tr>
  <tr> <td> T&gt;G </td> <td align="right"> 262510 </td> </tr>
   </table>

### INDELs
#### BRCA1
##### Count


```r
query <- "./sql/indel-count.sql"
where <- list("#_AND_"="AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, where))
```

```
# Count total number of indels
SELECT
  COUNT(length_alt) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype,
    length_alt,
    length_ref,
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      LENGTH(alternate_bases) AS length_alt,
      LENGTH(reference_bases) AS length_ref,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD IF EVERY(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND LENGTH(reference_bases) > 1
      OR LENGTH(alternate_bases) > 1
  )
  WHERE
    genotype != '-1,-1'
    AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype,
    length_ref,
    length_alt)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:33 2015 -->
<table border=1>
<tr> <th> count </th>  </tr>
  <tr> <td align="right">  44 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/mutation_counts_brca1.csv'
expectedResult = read.csv(file)
bcftoolsCount = expectedResult[expectedResult$type=='indels',]
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:33 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> indels </td> <td align="right">  41 </td> </tr>
   </table>

##### Length Distibution


```r
query <- "./sql/indel-length-distribution.sql"
where <- list("#_AND_"="AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(queryReplacements, where))
```

```
# Count the number of indels by length
SELECT
  length_alt,
  length_ref,
  COUNT(length_alt) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype,
    length_alt,
    length_ref,
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      LENGTH(alternate_bases) AS length_alt,
      LENGTH(reference_bases) AS length_ref,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD IF EVERY (call.genotype <= 0)
    HAVING
      num_alts = 1
      AND LENGTH(reference_bases) > 1
      OR LENGTH(alternate_bases) > 1
  )
  WHERE 
    genotype != '-1,-1'
    AND reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype,
    length_ref,
    length_alt)
GROUP EACH BY 
  length_alt,
  length_ref,
ORDER BY 
  length_ref,
  length_alt DESC
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:36 2015 -->
<table border=1>
<tr> <th> length_alt </th> <th> length_ref </th> <th> count </th>  </tr>
  <tr> <td align="right">   7 </td> <td align="right">   1 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   6 </td> <td align="right">   1 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   5 </td> <td align="right">   1 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   4 </td> <td align="right">   1 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">   3 </td> <td align="right">   1 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">   2 </td> <td align="right">   1 </td> <td align="right">  16 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   2 </td> <td align="right">  11 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   3 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   6 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   7 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   8 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">  13 </td> <td align="right">   1 </td> </tr>
   </table>

Compare to bcftools

```r
file = './data/bcfstats/indel_length_brca1.csv'
expectedResult = read.csv(file)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:36 2015 -->
<table border=1>
<tr> <th> length </th> <th> count </th>  </tr>
  <tr> <td align="right"> -12 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  -7 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  -6 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  -5 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  -3 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">  -2 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  -1 </td> <td align="right">  11 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">  16 </td> </tr>
  <tr> <td align="right">   2 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">   3 </td> <td align="right">   3 </td> </tr>
  <tr> <td align="right">   4 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   5 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   6 </td> <td align="right">   1 </td> </tr>
   </table>

#### Full Genome


```r
query <- "./sql/indel-count.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count total number of indels
SELECT
  COUNT(length_alt) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype,
    length_alt,
    length_ref,
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      LENGTH(alternate_bases) AS length_alt,
      LENGTH(reference_bases) AS length_ref,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD IF EVERY(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND LENGTH(reference_bases) > 1
      OR LENGTH(alternate_bases) > 1
  )
  WHERE
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype,
    length_ref,
    length_alt)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:39 2015 -->
<table border=1>
<tr> <th> count </th>  </tr>
  <tr> <td align="right"> 1522513 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/mutation_counts_brca1.csv'
expectedResult = read.csv(file)
bcftoolsCount = expectedResult[expectedResult$type=='indels',]
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:39 2015 -->
<table border=1>
<tr> <th> type </th> <th> count </th>  </tr>
  <tr> <td> indels </td> <td align="right">  41 </td> </tr>
   </table>

##### Length Distibution


```r
query <- "./sql/indel-length-distribution.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Count the number of indels by length
SELECT
  length_alt,
  length_ref,
  COUNT(length_alt) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype,
    length_alt,
    length_ref,
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      LENGTH(alternate_bases) AS length_alt,
      LENGTH(reference_bases) AS length_ref,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
    OMIT RECORD IF EVERY (call.genotype <= 0)
    HAVING
      num_alts = 1
      AND LENGTH(reference_bases) > 1
      OR LENGTH(alternate_bases) > 1
  )
  WHERE 
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype,
    length_ref,
    length_alt)
GROUP EACH BY 
  length_alt,
  length_ref,
ORDER BY 
  length_ref,
  length_alt DESC
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:42 2015 -->
<table border=1>
<tr> <th> length_alt </th> <th> length_ref </th> <th> count </th>  </tr>
  <tr> <td align="right">  29 </td> <td align="right">   1 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">  28 </td> <td align="right">   1 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right">  27 </td> <td align="right">   1 </td> <td align="right">   8 </td> </tr>
  <tr> <td align="right">  26 </td> <td align="right">   1 </td> <td align="right">  25 </td> </tr>
  <tr> <td align="right">  25 </td> <td align="right">   1 </td> <td align="right">  49 </td> </tr>
  <tr> <td align="right">  24 </td> <td align="right">   1 </td> <td align="right">  61 </td> </tr>
   </table>

Compare to bcftools stats


```r
file = './data/bcfstats/indel_length.csv'
expectedResult = read.csv(file)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Apr  3 09:45:42 2015 -->
<table border=1>
<tr> <th> length </th> <th> count </th>  </tr>
  <tr> <td align="right"> -43 </td> <td align="right">   2 </td> </tr>
  <tr> <td align="right"> -42 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right"> -41 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right"> -40 </td> <td align="right">   9 </td> </tr>
  <tr> <td align="right"> -39 </td> <td align="right">  10 </td> </tr>
  <tr> <td align="right"> -38 </td> <td align="right">  27 </td> </tr>
  <tr> <td align="right"> -37 </td> <td align="right">  15 </td> </tr>
  <tr> <td align="right"> -36 </td> <td align="right">  50 </td> </tr>
  <tr> <td align="right"> -35 </td> <td align="right">  29 </td> </tr>
  <tr> <td align="right"> -34 </td> <td align="right">  86 </td> </tr>
  <tr> <td align="right"> -33 </td> <td align="right">  37 </td> </tr>
  <tr> <td align="right"> -32 </td> <td align="right"> 165 </td> </tr>
  <tr> <td align="right"> -31 </td> <td align="right">  73 </td> </tr>
  <tr> <td align="right"> -30 </td> <td align="right"> 238 </td> </tr>
  <tr> <td align="right"> -29 </td> <td align="right"> 107 </td> </tr>
  <tr> <td align="right"> -28 </td> <td align="right"> 382 </td> </tr>
  <tr> <td align="right"> -27 </td> <td align="right"> 194 </td> </tr>
  <tr> <td align="right"> -26 </td> <td align="right"> 539 </td> </tr>
  <tr> <td align="right"> -25 </td> <td align="right"> 293 </td> </tr>
  <tr> <td align="right"> -24 </td> <td align="right"> 984 </td> </tr>
  <tr> <td align="right"> -23 </td> <td align="right"> 389 </td> </tr>
  <tr> <td align="right"> -22 </td> <td align="right"> 1012 </td> </tr>
  <tr> <td align="right"> -21 </td> <td align="right"> 585 </td> </tr>
  <tr> <td align="right"> -20 </td> <td align="right"> 1986 </td> </tr>
  <tr> <td align="right"> -19 </td> <td align="right"> 759 </td> </tr>
  <tr> <td align="right"> -18 </td> <td align="right"> 2119 </td> </tr>
  <tr> <td align="right"> -17 </td> <td align="right"> 986 </td> </tr>
  <tr> <td align="right"> -16 </td> <td align="right"> 3727 </td> </tr>
  <tr> <td align="right"> -15 </td> <td align="right"> 2118 </td> </tr>
  <tr> <td align="right"> -14 </td> <td align="right"> 3939 </td> </tr>
  <tr> <td align="right"> -13 </td> <td align="right"> 2118 </td> </tr>
  <tr> <td align="right"> -12 </td> <td align="right"> 8764 </td> </tr>
  <tr> <td align="right"> -11 </td> <td align="right"> 2741 </td> </tr>
  <tr> <td align="right"> -10 </td> <td align="right"> 8215 </td> </tr>
  <tr> <td align="right">  -9 </td> <td align="right"> 4426 </td> </tr>
  <tr> <td align="right">  -8 </td> <td align="right"> 16036 </td> </tr>
  <tr> <td align="right">  -7 </td> <td align="right"> 5321 </td> </tr>
  <tr> <td align="right">  -6 </td> <td align="right"> 19232 </td> </tr>
  <tr> <td align="right">  -5 </td> <td align="right"> 18382 </td> </tr>
  <tr> <td align="right">  -4 </td> <td align="right"> 71630 </td> </tr>
  <tr> <td align="right">  -3 </td> <td align="right"> 41525 </td> </tr>
  <tr> <td align="right">  -2 </td> <td align="right"> 109592 </td> </tr>
  <tr> <td align="right">  -1 </td> <td align="right"> 281414 </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right"> 304088 </td> </tr>
  <tr> <td align="right">   2 </td> <td align="right"> 109146 </td> </tr>
  <tr> <td align="right">   3 </td> <td align="right"> 33188 </td> </tr>
  <tr> <td align="right">   4 </td> <td align="right"> 61497 </td> </tr>
  <tr> <td align="right">   5 </td> <td align="right"> 14112 </td> </tr>
  <tr> <td align="right">   6 </td> <td align="right"> 15290 </td> </tr>
  <tr> <td align="right">   7 </td> <td align="right"> 3621 </td> </tr>
  <tr> <td align="right">   8 </td> <td align="right"> 11632 </td> </tr>
  <tr> <td align="right">   9 </td> <td align="right"> 2994 </td> </tr>
  <tr> <td align="right">  10 </td> <td align="right"> 4646 </td> </tr>
  <tr> <td align="right">  11 </td> <td align="right"> 1402 </td> </tr>
  <tr> <td align="right">  12 </td> <td align="right"> 4365 </td> </tr>
  <tr> <td align="right">  13 </td> <td align="right"> 978 </td> </tr>
  <tr> <td align="right">  14 </td> <td align="right"> 1470 </td> </tr>
  <tr> <td align="right">  15 </td> <td align="right"> 950 </td> </tr>
  <tr> <td align="right">  16 </td> <td align="right"> 1287 </td> </tr>
  <tr> <td align="right">  17 </td> <td align="right"> 469 </td> </tr>
  <tr> <td align="right">  18 </td> <td align="right"> 568 </td> </tr>
  <tr> <td align="right">  19 </td> <td align="right"> 272 </td> </tr>
  <tr> <td align="right">  20 </td> <td align="right"> 351 </td> </tr>
  <tr> <td align="right">  21 </td> <td align="right"> 177 </td> </tr>
  <tr> <td align="right">  22 </td> <td align="right"> 104 </td> </tr>
  <tr> <td align="right">  23 </td> <td align="right">  60 </td> </tr>
  <tr> <td align="right">  24 </td> <td align="right">  46 </td> </tr>
  <tr> <td align="right">  25 </td> <td align="right">  24 </td> </tr>
  <tr> <td align="right">  26 </td> <td align="right">   8 </td> </tr>
  <tr> <td align="right">  27 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right">  28 </td> <td align="right">   1 </td> </tr>
   </table>


