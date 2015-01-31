<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->
  
  <!-- Copyright 2015 Stanford University All rights reserved. -->
  
  Performing QC on AAA genomes using Google Genomics
================================================
  
  The following example makes use of 5 genomes from the MVP.

Setting Up and Describing the Data
----------------------------------
  




```r
# Setup for BigQuery access
require(bigrquery)
require(xtable)
require(RCurl)
require(dplyr)

project <- "gbsc-gcp-project-mvp"                   # put your projectID here

DisplayAndDispatchQuery <- function(queryUri, replacements=list()) {
  if(grepl("^https.*", queryUri)) {
    querySql <- getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    querySql <- readChar(queryUri, nchars=1e6)
  }
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }
  cat(querySql)
  query_exec(querySql, project)
}

table_replacement <- list("_THE_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_gvcfs",
                          "_THE_EXPANDED_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_vcfs")
```

Check Singletons
-----------------------------------


```r
limits = "WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30"
result <- DisplayAndDispatchQuery("./sql/private-variants-brca1.sql",
                                  replacements=c(table_replacement, "_LIMITS_" = limits))
```

```
# Private variants within BRCA1.
SELECT
  reference_name AS CHROM,
  start AS POS,
  GROUP_CONCAT(CASE WHEN cnt = 1 THEN 'S'
    WHEN cnt = 2 THEN 'D'
    ELSE STRING(cnt) END) AS SINGLETON_DOUBLETON,
  reference_bases AS REF,
  alternate_bases AS ALT,
  GROUP_CONCAT(call.call_set_name) AS INDV,
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
          FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_gvcfs], call)
        WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30
        OMIT
          call IF EVERY(call.genotype = -1)
          ),
        alternate_bases)
      )
  OMIT
    RECORD IF alternate_bases IS NULL
  HAVING
    cnt > 0
    )
GROUP EACH BY
  chrom,
  pos,
  ref,
  alt
HAVING
  num_samples_with_variant = 1
ORDER BY
  POS,
  INDV

Retrieving data:  3.6s
Retrieving data:  5.4s
Retrieving data:  7.1s
```
Number of rows returned by this query: 40418.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Jan 30 16:23:09 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALT </th> <th> INDV </th> <th> genotype </th> <th> num_samples_with_variant </th>  </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> D </td> <td> A </td> <td> C </td> <td> LP6005038-DNA_A03 </td> <td> "1,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050251 </td> <td> S </td> <td> A </td> <td> T </td> <td> LP6005038-DNA_A03 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050611 </td> <td> S </td> <td> C </td> <td> G </td> <td> LP6005038-DNA_B02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050821 </td> <td> S </td> <td> G </td> <td> A </td> <td> LP6005038-DNA_B02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050839 </td> <td> S </td> <td> C </td> <td> G </td> <td> LP6005038-DNA_A02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16051967 </td> <td> S </td> <td> C </td> <td> A </td> <td> LP6005038-DNA_B02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
   </table>

Check Individual Heterozygosity
-----------------------------------


```r
limits = "WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30"
result <- DisplayAndDispatchQuery("./sql/homozygous-variants.sql",
                                  replacements=c(table_replacement, "_LIMITS_"=limits))
```

```
# Individual Homozygosity
SELECT
  INDV,
  O_HOM,
  ROUND(E_HOM, 2) as E_HOM,
  N_SITES,
  ROUND((O_HOM - E_HOM) / (N_SITES - E_HOM), 5) AS F
FROM (
  SELECT
    call.call_set_name AS INDV,
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
      NTH(1,
        call.genotype) WITHIN call AS first_allele,
      NTH(2,
        call.genotype) WITHIN call AS second_allele,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype >= 0) WITHIN RECORD AS called_allele_count,
      IF((SUM(1 = call.genotype) > 0),
        SUM(call.genotype = 1)/SUM(call.genotype >= 0),
        -1)  WITHIN RECORD AS freq
    FROM
      FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_vcfs], call)
    WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30
    OMIT
      call IF SOME(call.genotype < 0)
      OR (2 > COUNT(call.genotype))
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    INDV
    )
ORDER BY
  INDV
```
Number of rows returned by this query: 5.

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Jan 30 16:23:14 2015 -->
<table border=1>
<tr> <th> INDV </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 63849 </td> <td align="right"> 444361.00 </td> <td align="right"> 97599 </td> <td align="right"> 1.10 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 64435 </td> <td align="right"> 444939.00 </td> <td align="right"> 97574 </td> <td align="right"> 1.10 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 63852 </td> <td align="right"> 445228.00 </td> <td align="right"> 97681 </td> <td align="right"> 1.10 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 64104 </td> <td align="right"> 449008.00 </td> <td align="right"> 97710 </td> <td align="right"> 1.10 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 64855 </td> <td align="right"> 442527.00 </td> <td align="right"> 97617 </td> <td align="right"> 1.09 </td> </tr>
   </table>

Cohort Level QC
===============

Check Hardy-Weinberg Equilibrium
-----------------------------------

```r
limits = "WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30"
result <- DisplayAndDispatchQuery("./sql/hardy-weinberg-brca1-expanded.sql",
                                  replacements=c(table_replacement, "_LIMITS_"=limits))
```

```
# The following query computes the Hardy-Weinberg equilibrium for BRCA1 variants.
SELECT
  CHR,
  POS,
  ref,
  alt,
  OBS_HOM1,
  OBS_HET,
  OBS_HOM2,
  E_HOM1,
  E_HET,
  E_HOM2,

  # Chi Squared Calculation
  # SUM(((Observed - Expected)^2) / Expected )
  ROUND((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2), 6)
  AS ChiSq,

  # Determine if Chi Sq value is significant
  IF((POW(OBS_HOM1 - E_HOM1, 2) / E_HOM1)
  + (POW(OBS_HET - E_HET, 2) / E_HET)
  + (POW(OBS_HOM2 - E_HOM2, 2) / E_HOM2)
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG

FROM (
    SELECT
      CHR,
      POS,
      ref,
      alt,
      OBS_HOM1,
      OBS_HET,
      OBS_HOM2,

      # Expected AA
      # p^2
      # ((COUNT(AA) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
      ROUND(POW((OBS_HOM1 + (OBS_HET/2)) /
        SAMPLE_COUNT, 2) * SAMPLE_COUNT, 2)
        AS E_HOM1,

      # Expected Aa
      # 2pq
      # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) *
      # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT)
      # * SAMPLE_COUNT
      ROUND(2 * ((OBS_HOM1 + (OBS_HET/2)) / SAMPLE_COUNT) *
        ((OBS_HOM2 + (OBS_HET/2)) / SAMPLE_COUNT)
        * SAMPLE_COUNT, 2)
        AS E_HET,

      # Expected aa
      # q^2
      # (COUNT(aa) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT
      ROUND(POW((OBS_HOM2 + (OBS_HET/2)) /
        SAMPLE_COUNT, 2) * SAMPLE_COUNT, 2)
        AS E_HOM2,

  FROM (
    SELECT
      reference_name AS CHR,
      start AS POS,
      reference_bases AS ref,
      alternate_bases AS alt,
      HOM_REF AS OBS_HOM1,
      HET AS OBS_HET,
      HOM_ALT AS OBS_HOM2,
      HOM_REF + HET + HOM_ALT AS SAMPLE_COUNT,
    FROM (
      SELECT
        reference_name,
        start,
        END,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
        SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
        SUM(SOME(0 = call.genotype)
          AND SOME(1 = call.genotype)) WITHIN call AS HET,
      FROM
        FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_vcfs], call)
      WHERE 
          reference_name = 'chr22'
          AND call.QUAL >= 30
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        )))
ORDER BY
  CHR,
  POS,
  ref

Running query:   RUNNING  2.7s
Running query:   RUNNING  3.4s
Running query:   RUNNING  4.0s
Running query:   RUNNING  4.7s
Running query:   RUNNING  5.3s
Running query:   RUNNING  5.9s
Running query:   RUNNING  6.6s
Running query:   RUNNING  7.3s
Running query:   RUNNING  7.9s

Retrieving data:  3.2s
Retrieving data:  5.6s
Retrieving data:  7.7s
Retrieving data: 11.0s
Retrieving data: 13.4s
Retrieving data: 15.9s
Retrieving data: 18.1s
Retrieving data: 20.9s
Retrieving data: 23.7s
```
Number of rows returned by this query: 100000.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Jan 30 16:23:51 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> A </td> <td> C </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> A </td> <td> C </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> A </td> <td> C </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> A </td> <td> C </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050035 </td> <td> A </td> <td> C </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 1.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr22 </td> <td align="right"> 16050158 </td> <td> C </td> <td> T </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
   </table>

===============

Check Transition-Transversion Ratio
-----------------------------------

```r
limits = "AND reference_name = 'chr22'
          AND call.QUAL >= 30"
result <- DisplayAndDispatchQuery("./sql/ti-tv-ratio.sql",
                                  replacements=c(table_replacement,"_LIMITS_"=limits))
```

```
# Compute the Ti/Tv ratio of the 1,000 Genomes dataset.
SELECT
  transitions,
  transversions,
  transitions/transversions AS titv,
  COUNT
FROM (
  SELECT
    SUM(IF(mutation IN ('A->G',
          'G->A',
          'C->T',
          'T->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transitions,
    SUM(IF(mutation IN ('A->C',
          'C->A',
          'G->T',
          'T->G',
          'A->T',
          'T->A',
          'C->G',
          'G->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transversions,
        COUNT(mutation) AS COUNT
  FROM (
    SELECT
      CONCAT(reference_bases,
        CONCAT(STRING('->'),
          alternate_bases)) AS mutation,
      COUNT(alternate_bases) AS num_snps,
    FROM
      FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.sample_gvcfs], call)
    WHERE
      LENGTH(alternate_bases) == 1
      AND LENGTH(reference_bases) == 1
      AND reference_name = 'chr22'
          AND call.QUAL >= 30
    GROUP BY
      mutation,
    ORDER BY
      mutation))
```
The result:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Jan 30 16:23:53 2015 -->
<table border=1>
<tr> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> COUNT </th>  </tr>
  <tr> <td align="right"> 174857 </td> <td align="right"> 75357 </td> <td align="right"> 2.32 </td> <td align="right">  12 </td> </tr>
   </table>

