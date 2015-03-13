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
  
# Comparison of Google Genomics tools vs Academic Tools
  






```r
tableReplacement <- list("_THE_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls",
                          "_THE_EXPANDED_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls")
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

constraints <- list("#_WHERE_"="WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")

googleRepository = "https://raw.githubusercontent.com/deflaux/codelabs/qc-codelab/R/PlatinumGenomes-QC/sql/"
```
The purpose of this codelab is to compare BigQuery queries to standard command line tools for genomic analysis.  

* [Singleton Rate](#singleton-rate)

## Singleton Rate

For each genome, count the number of variants shared by no other member of the cohort.  Too many private calls for a particular individual may indicate a problem.


```r
query <- paste(googleRepository, "private-variants-brca1.sql", sep = "")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=tableReplacement)
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
  COUNT(call.call_set_name) AS num_samples_with_variant
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
        WHERE
          reference_name = 'chr17'
          AND start BETWEEN 41196311
          AND 41277499
        OMIT call IF EVERY(call.genotype = -1)
          ),
        alternate_bases)
      )
  OMIT RECORD IF alternate_bases IS NULL
  HAVING
    cnt > 0
    )
GROUP BY
  chrom,
  pos,
  ref,
  alt
HAVING
  num_samples_with_variant = 1
ORDER BY
  POS,
  INDV
```
Number of rows returned by this query: 33.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:51:50 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALT </th> <th> INDV </th> <th> genotype </th> <th> num_samples_with_variant </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944 </td> <td> S </td> <td> T </td> <td> C </td> <td> LP6005038-DNA_B02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938 </td> <td> S </td> <td> A </td> <td> AT </td> <td> LP6005038-DNA_B02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td> S </td> <td> T </td> <td> C </td> <td> LP6005038-DNA_A02 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202633 </td> <td> S </td> <td> AG </td> <td> A </td> <td> LP6005038-DNA_A01 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204835 </td> <td> S </td> <td> ACACACT </td> <td> A </td> <td> LP6005038-DNA_A01 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208691 </td> <td> S </td> <td> G </td> <td> GA </td> <td> LP6005038-DNA_A03 </td> <td> "0,1" </td> <td align="right">   1 </td> </tr>
   </table>

Compare to [output](./data/brca1/singletons.singletons) from vcftools, which has 85 some of which are for 0/0 genotypes from reference matching blocks (see the [vcftools command line](./data/brca1/singletons.e9636222) used to create this file).


```r
expectedResult <- read.table("./data/brca1/singletons.singletons", header=TRUE)
# Convert to zero-based coordinates
expectedResult <- mutate(expectedResult, POS = POS - 1)
# Clean colnames to match
colnames(expectedResult) <- gsub('\\.+', '_', colnames(expectedResult))
```

How many singletons do the two results have in common?

```r
nrow(inner_join(result, expectedResult))
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```
## [1] 54
```

Which singletons were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
print(xtable(onlyBQ), type="html", include.rownames=F)
```

```
## Warning in matrix(align.tmp[(2 - pos):(ncol(x) + 1)], nrow = nrow(x), ncol
## = ncol(x) + : data length exceeds size of matrix
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:51:50 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> REF </th> <th> ALT </th> <th> INDV </th> <th> genotype </th> <th> num_samples_with_variant </th>  </tr>
  </table>

Which singletons were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
print(xtable(onlyVcftools), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:51:50 2015 -->
<table border=1>
<tr> <th> CHROM </th> <th> POS </th> <th> SINGLETON_DOUBLETON </th> <th> ALLELE </th> <th> INDV </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271293.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267518.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256102.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256101.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256100.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256086.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256085.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256084.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256083.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256082.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256081.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256080.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256079.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256078.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256077.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256076.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256075.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252696.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B01 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252647.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252634.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249363.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242077.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242076.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242075.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230105.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229760.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226740.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226739.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226738.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226737.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226736.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214210.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214209.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202634.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
   </table>

Retrieving the gVCF data for the singletons identified only by vcftools:

```r
having <- paste("start = ", onlyVcftools$POS,
                sep="", collapse=" OR ")
query <- paste(googleRepository, "examine-data.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement,
                                                 "_HAVING_"=having))
```

```
# Examine the data for particular calls.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
FROM
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
WHERE
  reference_name = 'chr17'
HAVING
  start = 41271293 OR start = 41267518 OR start = 41256102 OR start = 41256101 OR start = 41256100 OR start = 41256086 OR start = 41256085 OR start = 41256084 OR start = 41256083 OR start = 41256082 OR start = 41256081 OR start = 41256080 OR start = 41256079 OR start = 41256078 OR start = 41256077 OR start = 41256076 OR start = 41256075 OR start = 41252696 OR start = 41252647 OR start = 41252634 OR start = 41249363 OR start = 41242077 OR start = 41242076 OR start = 41242075 OR start = 41230105 OR start = 41229760 OR start = 41226740 OR start = 41226739 OR start = 41226738 OR start = 41226737 OR start = 41226736 OR start = 41214210 OR start = 41214209 OR start = 41202634
ORDER BY
  start,
  end,
  call.call_set_name
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:51:54 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> genotype </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th>  </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41252696 </td> <td align="right"> 41252755 </td> <td> T </td> <td>  </td> <td align="right"> 31.21 </td> <td>  </td> </tr>
   </table>

Most of these positions correspond to sites where 4 of the 5 genomes had a no-call, this is not neccessarily a singleton.

It appears that they correspond either to:
* a reference-matching block, so not actually a singleton and just perhaps violating an assumption in the vcftools code
* or a non-singleon variant, perhaps due to a problem in converting the gVCF data to all-positions VCF via gvcftools?


## Homozygosity Rate and Inbreeding Coefficient

For each genome, compare the expected and observed rates of homozygosity.


```r
query <- paste(googleRepository, "homozygous-variants.sql", sep="") 
result <- DisplayAndDispatchQuery("./sql/homozygous-variants.sql",
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
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
    WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
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
Number of rows returned by this query: 5.

Displaying the first few results:
  <!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
  <!-- Fri Mar 13 10:51:56 2015 -->
  <table border=1>
  <tr> <th> call_call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
    <tr> <td> LP6005038-DNA_A01 </td> <td align="right">  17 </td> <td align="right"> 80.31 </td> <td align="right"> 137 </td> <td align="right"> -1.12 </td> </tr>
    <tr> <td> LP6005038-DNA_A02 </td> <td align="right">  16 </td> <td align="right"> 75.11 </td> <td align="right"> 137 </td> <td align="right"> -0.96 </td> </tr>
    <tr> <td> LP6005038-DNA_A03 </td> <td align="right">  17 </td> <td align="right"> 74.48 </td> <td align="right"> 136 </td> <td align="right"> -0.93 </td> </tr>
    <tr> <td> LP6005038-DNA_B01 </td> <td align="right">  20 </td> <td align="right"> 86.78 </td> <td align="right"> 138 </td> <td align="right"> -1.30 </td> </tr>
    <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 141 </td> <td align="right"> 115.11 </td> <td align="right"> 144 </td> <td align="right"> 0.90 </td> </tr>
     </table>



```r
expectedResult <- read.table("./data/brca1/het_calls.het", header=TRUE)
# Clean colnames to match
colnames(expectedResult) <- gsub('\\.+$', '', colnames(expectedResult))
colnames(expectedResult) <- gsub('\\.+', '_', colnames(expectedResult))
```


```r
n = names(result)
n[1] = "INDV"
names(result) = n
joinedResult <- inner_join(expectedResult, result, by=c("INDV"))
print(xtable(joinedResult[,order(colnames(joinedResult))]), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:51:57 2015 -->
<table border=1>
<tr> <th> E_HOM.x </th> <th> E_HOM.y </th> <th> F.x </th> <th> F.y </th> <th> INDV </th> <th> N_SITES.x </th> <th> N_SITES.y </th> <th> O_HOM.x </th> <th> O_HOM.y </th>  </tr>
  <tr> <td align="right"> 65.50 </td> <td align="right"> 80.31 </td> <td align="right"> -0.85 </td> <td align="right"> -1.12 </td> <td> LP6005038-DNA_A01 </td> <td align="right"> 143 </td> <td align="right"> 137 </td> <td align="right">   0 </td> <td align="right">  17 </td> </tr>
  <tr> <td align="right"> 66.40 </td> <td align="right"> 75.11 </td> <td align="right"> -0.70 </td> <td align="right"> -0.96 </td> <td> LP6005038-DNA_A02 </td> <td align="right"> 150 </td> <td align="right"> 137 </td> <td align="right">   8 </td> <td align="right">  16 </td> </tr>
  <tr> <td align="right"> 67.60 </td> <td align="right"> 74.48 </td> <td align="right"> -0.72 </td> <td align="right"> -0.93 </td> <td> LP6005038-DNA_A03 </td> <td align="right"> 149 </td> <td align="right"> 136 </td> <td align="right">   9 </td> <td align="right">  17 </td> </tr>
  <tr> <td align="right"> 64.00 </td> <td align="right"> 86.78 </td> <td align="right"> -0.87 </td> <td align="right"> -1.30 </td> <td> LP6005038-DNA_B01 </td> <td align="right"> 135 </td> <td align="right"> 138 </td> <td align="right">   2 </td> <td align="right">  20 </td> </tr>
  <tr> <td align="right"> 61.80 </td> <td align="right"> 115.11 </td> <td align="right"> 0.92 </td> <td align="right"> 0.90 </td> <td> LP6005038-DNA_B02 </td> <td align="right"> 133 </td> <td align="right"> 144 </td> <td align="right"> 127 </td> <td align="right"> 141 </td> </tr>
   </table>

The logic in the query looks similar to vcftools [output_het method](http://sourceforge.net/p/vcftools/code/HEAD/tree/trunk/cpp/variant_file_output.cpp#l165) but there is clearly a difference.  TODO: investigate the difference further.


Check Hardy-Weinberg Equilibrium
-----------------------------------

```r
query <- paste(googleRepository, "hardy-weinberg-brca1.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=tableReplacement)
```

```
# The following query computes the Hardy-Weinberg equilibrium for BRCA1 SNPs.
SELECT
  calcs.CHR AS CHR,
  calcs.POS AS POS,
  calcs.ref AS ref,
  calcs.alt AS alt,
  calcs.OBS_HOM1 AS OBS_HOM1,
  calcs.OBS_HET AS OBS_HET,
  calcs.OBS_HOM2 AS OBS_HOM2,
  calcs.EXP_HOM1 AS EXP_HOM1,
  calcs.EXP_HET AS EXP_HET,
  calcs.EXP_HOM2 AS EXP_HOM2,

  # Chi Squared Calculation
  # SUM(((Observed - Expected)^2) / Expected )
  ROUND((POW(calcs.OBS_HOM1 - calcs.EXP_HOM1, 2) / calcs.EXP_HOM1)
  + (POW(calcs.OBS_HET - calcs.EXP_HET, 2) / calcs.EXP_HET)
  + (POW(calcs.OBS_HOM2 - calcs.EXP_HOM2, 2) / calcs.EXP_HOM2), 3)
  AS CHI_SQ,

  # Determine if Chi Sq value is significant
  IF((POW(calcs.OBS_HOM1 - calcs.EXP_HOM1, 2) / calcs.EXP_HOM1)
  + (POW(calcs.OBS_HET - calcs.EXP_HET, 2) / calcs.EXP_HET)
  + (POW(calcs.OBS_HOM2 - calcs.EXP_HOM2, 2) / calcs.EXP_HOM2)
  > 5.991, "TRUE", "FALSE") AS PVALUE_SIG

FROM (
    SELECT
      vals.CHR AS CHR,
      vals.POS AS POS,
      vals.ref AS ref,
      vals.alt AS alt,
      vals.OBS_HOM1 AS OBS_HOM1,
      vals.OBS_HET AS OBS_HET,
      vals.OBS_HOM2 AS OBS_HOM2,

      # Expected AA
      # p^2
      # ((COUNT(AA) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2) * SAMPLE_COUNT
      ROUND(POW((vals.OBS_HOM1 + (vals.OBS_HET/2)) /
        vals.SAMPLE_COUNT, 2) * vals.SAMPLE_COUNT, 2)
        AS EXP_HOM1,

      # Expected Aa
      # 2pq
      # 2 * (COUNT(AA) + (COUNT(Aa)/2) / SAMPLE_COUNT) *
      # (COUNT(aa) + (COUNT(Aa)/2) / SAMPLE_COUNT)
      # * SAMPLE_COUNT
      ROUND(2 * ((vals.OBS_HOM1 + (vals.OBS_HET/2)) / vals.SAMPLE_COUNT) *
        ((vals.OBS_HOM2 + (vals.OBS_HET/2)) / vals.SAMPLE_COUNT)
        * vals.SAMPLE_COUNT, 2)
        AS EXP_HET,

      # Expected aa
      # q^2
      # (COUNT(aa) + (COUNT(Aa)/2) /
      #  SAMPLE_COUNT) ^ 2 * SAMPLE_COUNT
      ROUND(POW((vals.OBS_HOM2 + (vals.OBS_HET/2)) /
        vals.SAMPLE_COUNT, 2) * vals.SAMPLE_COUNT, 2)
        AS EXP_HOM2,

    FROM (
        SELECT
          vars.reference_name AS CHR,
          vars.start AS POS,
          reference_bases AS ref,
          alternate_bases AS alt,
          SUM(refs.HOM_REF) + vars.HOM_REF AS OBS_HOM1,
          vars.HET AS OBS_HET,
          vars.HOM_ALT AS OBS_HOM2,
          SUM(refs.HOM_REF) + vars.HOM_REF + vars.HET + vars.HOM_ALT AS SAMPLE_COUNT,

        FROM (
              # Constrain the left hand side of the _join to reference-matching blocks.
            SELECT
              reference_name,
              start,
              END,
              SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
            FROM
              [genomics-public-data:platinum_genomes.variants]
            WHERE
              reference_name = 'chr17'
            OMIT RECORD IF EVERY(alternate_bases IS NOT NULL)
              ) AS refs
          JOIN (
              SELECT
                reference_name,
                start,
                END,
                reference_bases,
                GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
                COUNT(alternate_bases) WITHIN RECORD AS num_alts,
                SUM(EVERY(0 = call.genotype)) WITHIN call AS HOM_REF,
                SUM(EVERY(1 = call.genotype)) WITHIN call AS HOM_ALT,
                SUM(SOME(0 = call.genotype) AND SOME(1 = call.genotype)) WITHIN call AS HET,

              FROM
                [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
              WHERE
                reference_name = 'chr17'
                AND start BETWEEN 41196311
                AND 41277499
            #  OMIT call IF 2 != COUNT(call.genotype)
              HAVING
                # Skip ref-matching blocks, 1/2 genotypes, and non-SNP variants
                num_alts = 1
                AND reference_bases IN ('A','C','G','T')
                AND alternate_bases IN ('A','C','G','T')
                ) AS vars
              # The _join criteria _is complicated since we are trying to see if a variant
              # overlaps a reference-matching interval.
            ON
              vars.reference_name = refs.reference_name
            WHERE
              refs.start <= vars.start
              AND refs.END >= vars.start+1
            GROUP BY
              CHR,
              POS,
              ref,
              alt,
              vars.HOM_REF,
              OBS_HET,
              OBS_HOM2,
              vars.HET,
              vars.HOM_ALT
            ORDER BY
              CHR,
              POS,
              ref,
              alt ) AS vals ) AS calcs
```
Number of rows returned by this query: 144.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:02 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> EXP_HOM1 </th> <th> EXP_HET </th> <th> EXP_HOM2 </th> <th> CHI_SQ </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198620 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41199912 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200108 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
   </table>

Compare to [brca1.hwe](./data/brca1/hardy.hwe) (see the [vcftools command line](./data/hwe/brca1.log) used to create this file).


```r
require(dplyr)
df <- read.table("./data/brca1/hardy.hwe", header=TRUE)
obsSplitCol <- "OBS.HOM1.HET.HOM2."
obsTemp <- read.table(text=as.character(df[, obsSplitCol]), sep = "/")
names(obsTemp) <- c("OBS_HOM1", "OBS_HET", "OBS_HOM2")
eSplitCol <- "E.HOM1.HET.HOM2."
eTemp <- read.table(text=as.character(df[, eSplitCol]), sep = "/")
names(eTemp) <- c("E_HOM1", "E_HET", "E_HOM2")
expectedResult <- cbind(cbind(df[setdiff(names(df), c(obsSplitCol,eSplitCol))], obsTemp), eTemp)
# Convert to zero-based coordinates
expectedResult <- mutate(expectedResult, POS = POS - 1)
```

How many results do the two results have in common?

```r
nrow(inner_join(result, expectedResult, by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2")))
```

```
## [1] 0
```

Which results were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult, , by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
print(xtable(arrange(onlyBQ, CHR, POS)), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:02 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> EXP_HOM1 </th> <th> EXP_HET </th> <th> EXP_HOM2 </th> <th> CHI_SQ </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198620 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41199912 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200108 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 10.02 </td> <td align="right"> 0.95 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201363 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201701 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202687 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203324 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203590 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204376 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204389 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204855 </td> <td> G </td> <td> A </td> <td align="right">   8 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 8.33 </td> <td align="right"> 3.33 </td> <td align="right"> 0.33 </td> <td align="right"> 0.48 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205771 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205940 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41206055 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41209577 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210395 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211652 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212168 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212337 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212546 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212804 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213625 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213659 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213747 </td> <td> T </td> <td> C </td> <td align="right">   8 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 8.33 </td> <td align="right"> 3.33 </td> <td align="right"> 0.33 </td> <td align="right"> 0.48 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213759 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213892 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213995 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41215824 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216204 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216205 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216932 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217873 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218332 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218571 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219340 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219559 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219779 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219803 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220222 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220287 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220771 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222098 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222461 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222722 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223093 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224832 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225780 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225782 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225838 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226600 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226674 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229385 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229811 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229856 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229907 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230227 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230335 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230375 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230523 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230536 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230989 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231220 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231515 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231901 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232697 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41234469 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41235798 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41237952 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239471 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239490 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239627 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239698 </td> <td> G </td> <td> A </td> <td align="right">  15 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 15.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239979 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41240276 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241389 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241502 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242284 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243189 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243999 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244434 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244935 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245236 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245465 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245470 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247054 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247603 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247614 </td> <td> A </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247623 </td> <td> A </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248163 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248483 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248903 </td> <td> G </td> <td> A </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249093 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250922 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251494 </td> <td> C </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251645 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251930 </td> <td> G </td> <td> A </td> <td align="right">   7 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 7.03 </td> <td align="right"> 0.94 </td> <td align="right"> 0.03 </td> <td align="right"> 0.03 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252574 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252611 </td> <td> T </td> <td> A </td> <td align="right">  11 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 11.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252612 </td> <td> T </td> <td> C </td> <td align="right">  11 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 11.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252614 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 17.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254173 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254404 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254485 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255101 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255110 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257133 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257457 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258042 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258152 </td> <td> G </td> <td> A </td> <td align="right">  13 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 13.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258790 </td> <td> T </td> <td> G </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259048 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260807 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261232 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263043 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263565 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264145 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264363 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41265775 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267049 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 10.02 </td> <td align="right"> 0.95 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270228 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270276 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270462 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270665 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273347 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273378 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273536 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274777 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275150 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275644 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276246 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276347 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41277186 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
   </table>

Note vcftools appears to skip variants with single allele genotypes:
```
zgrep 41242078 platinum_genomes_brca1_expanded_merged.vcf.gz 
chr17  41242078  .  G  A  143	LowGQX;TruthSensitivityTranche99.90to100.00;LowQD;SiteConflict	BLOCKAVG_min30p3a;MQ=57;MQ0=0;BaseQRankSum=0.781;Dels=0.3;FS=1.561;HRun=11;HaplotypeScore=77.7361;MQRankSum=0.093;QD=2.01;ReadPosRankSum=-2.871;SB=-45.67;VQSLOD=-1.8762;culprit=QD;set=FilteredInAll;DP=425;AF=0.5;AN=25;AC=1	GT:DP:GQX:MQ:AD:GQ:PL:VF	0/0:57:99:59:.:.:.:.	0:27:25:57:26:25.38:.:.	0/0:51:99:57:.:.:.:.	0/1:50:99:59:42,8:99:173,0,1238:0.16	0/0:46:99:59:.:.:.:.	0/0:44:99:60:.:.:.:.	.:46:.:59:40,6:.:.:.	0/0:40:85:59:.:.:.:.	0/0:40:85:59:.:.:.:.	0/0:63:99:58:.:.:.:.	0:42:2:58:37:1.58:.:.	0:33:0:57:29:0.03:.:.	.:44:.:58:31,12:.:.:.	0/0:44:90:58:.:.:.:.	0/0:40:87:58:.:.:.:.	.:44:.:57:39,5:.:.:.	0/0:55:99:59:.:.:.:.
```

Which results were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result, , by=c("CHR", "POS", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
print(xtable(arrange(onlyVcftools, CHR, POS)), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:02 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ChiSq </th> <th> P </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197938.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198620.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41199912.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200108.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536.00 </td> <td align="right"> 0.12 </td> <td align="right"> 1.00 </td> <td align="right">   2 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 2.08 </td> <td align="right"> 0.83 </td> <td align="right"> 0.08 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200549.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200703.00 </td> <td align="right"> 3.00 </td> <td align="right"> 0.40 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right"> 0.75 </td> <td align="right"> 1.50 </td> <td align="right"> 0.75 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201363.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201701.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202631.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202633.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202687.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203324.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203590.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204376.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204389.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204835.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204855.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205771.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205940.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41206055.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41208691.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41209577.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210395.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211652.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212168.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212337.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212546.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212804.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213625.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213659.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213747.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213759.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213892.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213995.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214208.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41215824.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216204.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216205.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216932.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217873.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218332.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218571.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219340.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219559.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219779.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219803.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220222.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220287.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220569.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220771.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222098.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222461.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222722.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223093.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224832.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224891.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225653.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225780.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225782.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225838.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226600.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226674.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226735.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41227082.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229351.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229385.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229759.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772.00 </td> <td align="right"> 0.75 </td> <td align="right"> 1.00 </td> <td align="right">   1 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 1.33 </td> <td align="right"> 1.33 </td> <td align="right"> 0.33 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229776.00 </td> <td align="right">  </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   5 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 5.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229811.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229856.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229907.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230104.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230227.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230335.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230375.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230523.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230536.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230989.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231220.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231515.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231901.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232343.00 </td> <td align="right">  </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   5 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 5.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232697.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41234469.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41235798.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41237952.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239471.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239490.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239627.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239914.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239979.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41240276.00 </td> <td align="right"> 0.31 </td> <td align="right"> 1.00 </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241389.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241502.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241567.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242074.00 </td> <td align="right"> 3.00 </td> <td align="right"> 0.40 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right"> 0.75 </td> <td align="right"> 1.50 </td> <td align="right"> 0.75 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242284.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243189.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243999.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244434.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244935.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245236.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245465.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245470.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247054.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247121.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247603.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247614.00 </td> <td align="right"> 0.06 </td> <td align="right"> 1.00 </td> <td align="right">   4 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 4.05 </td> <td align="right"> 0.90 </td> <td align="right"> 0.05 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247623.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248163.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248483.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248903.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249093.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249362.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250046.00 </td> <td align="right">  </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250677.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250922.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251494.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251645.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251930.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252230.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252574.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252590.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252611.00 </td> <td align="right"> 0.22 </td> <td align="right"> 1.00 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 1.12 </td> <td align="right"> 0.75 </td> <td align="right"> 0.12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252612.00 </td> <td align="right"> 0.22 </td> <td align="right"> 1.00 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 1.12 </td> <td align="right"> 0.75 </td> <td align="right"> 0.12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252694.00 </td> <td align="right"> 0.12 </td> <td align="right"> 1.00 </td> <td align="right">   2 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 2.08 </td> <td align="right"> 0.83 </td> <td align="right"> 0.08 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252705.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254173.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254404.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254485.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254964.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255101.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255110.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256074.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41256099.00 </td> <td align="right"> 3.00 </td> <td align="right"> 0.40 </td> <td align="right">   0 </td> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right"> 0.75 </td> <td align="right"> 1.50 </td> <td align="right"> 0.75 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257133.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257457.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258042.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258790.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259048.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260351.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260807.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261232.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263043.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263116.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263565.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264145.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264363.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264472.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41265775.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267049.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267517.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205.00 </td> <td align="right"> 0.22 </td> <td align="right"> 1.00 </td> <td align="right">   1 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 1.12 </td> <td align="right"> 0.75 </td> <td align="right"> 0.12 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270228.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270276.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270462.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270665.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270777.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271292.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273347.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273378.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273536.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274777.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905.00 </td> <td align="right"> 0.75 </td> <td align="right"> 1.00 </td> <td align="right">   1 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 1.33 </td> <td align="right"> 1.33 </td> <td align="right"> 0.33 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275080.00 </td> <td align="right"> 4.00 </td> <td align="right"> 0.31 </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275150.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275644.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276246.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276347.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41277186.00 </td> <td align="right"> 2.22 </td> <td align="right"> 0.43 </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> </tr>
   </table>

Retrieving the gVCF data for the results identified only by vcftools:

```r
having <- paste("start <= ", onlyVcftools$POS,
                "AND",
                "end >= ", onlyVcftools$POS+1)
query <- paste(googleRepository, "examine-data.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement,
                                                 "_HAVING_"=having))
```

```
# Examine the data for particular calls.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
FROM
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
WHERE
  reference_name = 'chr17'
HAVING
  _HAVING_
ORDER BY
  start,
  end,
  call.call_set_name
```

```
Error: Field '_HAVING_' not found.

query invalid. Field '_HAVING_' not found.
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:03 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> ref </th> <th> alt </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> EXP_HOM1 </th> <th> EXP_HET </th> <th> EXP_HOM2 </th> <th> CHI_SQ </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196407 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41197273 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41198620 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41199912 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200108 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 10.02 </td> <td align="right"> 0.95 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41200536 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201363 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201701 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41202687 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203324 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41203590 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204376 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204389 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41204855 </td> <td> G </td> <td> A </td> <td align="right">   8 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 8.33 </td> <td align="right"> 3.33 </td> <td align="right"> 0.33 </td> <td align="right"> 0.48 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205771 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41205940 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41206055 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41209577 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41210395 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211652 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212168 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212337 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212546 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41212804 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213625 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213659 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213747 </td> <td> T </td> <td> C </td> <td align="right">   8 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 8.33 </td> <td align="right"> 3.33 </td> <td align="right"> 0.33 </td> <td align="right"> 0.48 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213759 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213892 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213995 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41215824 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216204 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216205 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216932 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41217873 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218332 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41218571 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219340 </td> <td> G </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219559 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219779 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41219803 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220222 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220287 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220771 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222098 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222461 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222722 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41223093 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41224832 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225780 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225782 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41225838 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226600 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226674 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229385 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229772 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229811 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229856 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229907 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230227 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230335 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230375 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230523 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230536 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230989 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231220 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231515 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41231901 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41232697 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41234469 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41235798 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41237952 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239471 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239490 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239627 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239698 </td> <td> G </td> <td> A </td> <td align="right">  15 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 15.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41239979 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41240276 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241389 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41241502 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242284 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243189 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41243999 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244434 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41244935 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245236 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245465 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41245470 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 17.05 </td> <td align="right"> 1.89 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247054 </td> <td> C </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247603 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247614 </td> <td> A </td> <td> T </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41247623 </td> <td> A </td> <td> C </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248163 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248483 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41248903 </td> <td> G </td> <td> A </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41249093 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41250922 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251494 </td> <td> C </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251645 </td> <td> T </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41251930 </td> <td> G </td> <td> A </td> <td align="right">   7 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 7.03 </td> <td align="right"> 0.94 </td> <td align="right"> 0.03 </td> <td align="right"> 0.03 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252574 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252611 </td> <td> T </td> <td> A </td> <td align="right">  11 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 11.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252612 </td> <td> T </td> <td> C </td> <td align="right">  11 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 11.02 </td> <td align="right"> 0.96 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41252614 </td> <td> T </td> <td> C </td> <td align="right">  17 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 17.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254173 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254404 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41254485 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255101 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41255110 </td> <td> A </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257133 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41257457 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258042 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258152 </td> <td> G </td> <td> A </td> <td align="right">  13 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 13.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41258790 </td> <td> T </td> <td> G </td> <td align="right">  17 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 17.01 </td> <td align="right"> 0.97 </td> <td align="right"> 0.01 </td> <td align="right"> 0.01 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41259048 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41260807 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41261232 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263043 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41263565 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264145 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41264363 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41265775 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41267049 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 10.02 </td> <td align="right"> 0.95 </td> <td align="right"> 0.02 </td> <td align="right"> 0.02 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41268205 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270228 </td> <td> T </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270276 </td> <td> C </td> <td> T </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270462 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41270665 </td> <td> C </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273347 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273378 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273536 </td> <td> A </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274777 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 10.08 </td> <td align="right"> 1.83 </td> <td align="right"> 0.08 </td> <td align="right"> 0.10 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41274905 </td> <td> G </td> <td> A </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 10.00 </td> <td align="right"> 0.00 </td> <td align="right"> 0.00 </td> <td align="right">  </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275150 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41275644 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276246 </td> <td> A </td> <td> G </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41276347 </td> <td> T </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41277186 </td> <td> G </td> <td> C </td> <td align="right">  10 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 10.29 </td> <td align="right"> 3.43 </td> <td align="right"> 0.29 </td> <td align="right"> 0.39 </td> <td> FALSE </td> </tr>
   </table>

It appears that with BigQuery we are computing HWE for all the same variants as vcftools and the expected and Chi-Squared values are only slightly different.

See also: the [gVCF version of this query](./sql/hardy-weinberg-brca1.sql), which is close but only works for SNPs and needs a RIGHT OUTER JOIN to compute values for variants for which all the samples have the variant.

Ti/Tv by Alternate Allele Counts
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Variant-Level-QC.md#titv-by-alternate-allele-counts

```r
query <- paste(googleRepository, "ti-tv-by-alternate-allele-count.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
# Compute the Ti/Tv ratio for variants binned by alternate allele count.
SELECT
  transitions,
  transversions,
  transitions/transversions AS titv,
  alternate_allele_count
FROM (
  SELECT
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    alternate_allele_count
  FROM (
    SELECT
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      SUM(call.genotype = 1) WITHIN RECORD AS alternate_allele_count,
    FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_no_calls]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
  GROUP BY
    alternate_allele_count)
ORDER BY
  alternate_allele_count DESC
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:07 2015 -->
<table border=1>
<tr> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> alternate_allele_count </th>  </tr>
  <tr> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 0.00 </td> <td align="right">  10 </td> </tr>
  <tr> <td align="right">  81 </td> <td align="right">  29 </td> <td align="right"> 2.79 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right">  11 </td> <td align="right">   3 </td> <td align="right"> 3.67 </td> <td align="right">   2 </td> </tr>
  <tr> <td align="right">   8 </td> <td align="right">   5 </td> <td align="right"> 1.60 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   6 </td> <td align="right">   1 </td> <td align="right"> 6.00 </td> <td align="right">   0 </td> </tr>
   </table>
    
Compare to [output](./data/brca1/TsTv-by-count.TsTv.count) from vcftools.


```r
expectedResult <- read.table("./data/brca1/TsTv-by-count.TsTv.count", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:07 2015 -->
<table border=1>
<tr> <th> ALT_ALLELE_COUNT </th> <th> N_Ts </th> <th> N_Tv </th> <th> Ts.Tv </th>  </tr>
  <tr> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   1 </td> <td align="right">   8 </td> <td align="right">   5 </td> <td align="right"> 1.60 </td> </tr>
  <tr> <td align="right">   2 </td> <td align="right">  11 </td> <td align="right">   3 </td> <td align="right"> 3.67 </td> </tr>
  <tr> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   4 </td> <td align="right">  81 </td> <td align="right">  29 </td> <td align="right"> 2.79 </td> </tr>
  <tr> <td align="right">   5 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   6 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   7 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   8 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
  <tr> <td align="right">   9 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right">  </td> </tr>
   </table>




Missingness sample level
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Sample-Level-QC.md#missingness-rate

```r
query <- paste(googleRepository, "sample-level-missingness.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
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
    WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
    )
  GROUP BY
    call.call_set_name)
ORDER BY
  call.call_set_name
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:10 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 185 </td> <td align="right"> 81598 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right">  20 </td> <td align="right"> 81631 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 231 </td> <td align="right"> 81605 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 123 </td> <td align="right"> 81589 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right">  23 </td> <td align="right"> 82247 </td> <td align="right"> 0.00 </td> </tr>
   </table>
    
Compare to [output](./data/brca1/missingness.imiss) from vcftools.


```r
expectedResult <- read.table("./data/brca1/missingness.imiss", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:10 2015 -->
<table border=1>
<tr> <th> INDV </th> <th> N_DATA </th> <th> N_GENOTYPES_FILTERED </th> <th> N_MISS </th> <th> F_MISS </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 274 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 107 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 311 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 205 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right">  73 </td> <td align="right"> 0.00 </td> </tr>
   </table>

Missingness variant level
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Variant-Level-QC.md#missingness-rate

```r
query <- paste(googleRepository, "variant-level-missingness.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
# Compute the ratio no-calls for each variant.
SELECT
  reference_name,
  start,
  END,
  reference_bases,
  alternate_bases,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS missingness_rate
FROM (
  SELECT
    reference_name,
    start,
    END,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    SUM(call.genotype == -1) WITHIN RECORD AS no_calls,
    COUNT(call.genotype) WITHIN RECORD AS all_calls,
  FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  )
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:13 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> END </th> <th> reference_bases </th> <th> alternate_bases </th> <th> no_calls </th> <th> all_calls </th> <th> missingness_rate </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41235798 </td> <td align="right"> 41235799 </td> <td> G </td> <td> A </td> <td align="right">   0 </td> <td align="right">  10 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41271292 </td> <td align="right"> 41271294 </td> <td> GA </td> <td> G </td> <td align="right">   0 </td> <td align="right">   8 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41237952 </td> <td align="right"> 41237953 </td> <td> G </td> <td> A </td> <td align="right">   0 </td> <td align="right">  10 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273094 </td> <td align="right"> 41273095 </td> <td> G </td> <td> A </td> <td align="right">   0 </td> <td align="right">  10 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273347 </td> <td align="right"> 41273348 </td> <td> T </td> <td> C </td> <td align="right">   0 </td> <td align="right">  10 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41273378 </td> <td align="right"> 41273379 </td> <td> G </td> <td> C </td> <td align="right">   0 </td> <td align="right">  10 </td> <td align="right"> 0.00 </td> </tr>
   </table>
    
Compare to [output](./data/brca1/missing_sites.lmiss) from vcftools.


```r
expectedResult <- read.table("./data/brca1/missing_sites.lmiss", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:13 2015 -->
<table border=1>
<tr> <th> CHR </th> <th> POS </th> <th> N_DATA </th> <th> N_GENOTYPE_FILTERED </th> <th> N_MISS </th> <th> F_MISS </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196312 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196314 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196315 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196316 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196317 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
   </table>


Sex Inference
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Sample-Level-QC.md#sex-inference

```r
query <- paste(googleRepository, "gender-check.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
# Compute the the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the gender phenotype value is
# correct for each individual.
SELECT
  call.call_set_name,
  ROUND(SUM(het_RA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_het_alt_in_snvs,
  ROUND(SUM(hom_AA)/(SUM(hom_AA) + SUM(het_RA)), 3) AS perct_hom_alt_in_snvs,
  SUM(hom_AA) AS hom_AA_count,
  SUM(het_RA) AS het_RA_count,
  SUM(hom_RR) AS hom_RR_count,
FROM (
  SELECT
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    call.call_set_name,
    SOME(call.genotype = 0) AND NOT SOME(call.genotype > 0) WITHIN call AS hom_RR,
    SOME(call.genotype > 0) AND NOT SOME(call.genotype = 0) WITHIN call AS hom_AA,
    SOME(call.genotype > 0) AND SOME(call.genotype = 0) WITHIN call AS het_RA
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
  call.call_set_name
ORDER BY
  call.call_set_name
```

Number of rows returned by this query: 5.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri Mar 13 10:52:16 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> hom_AA_count </th> <th> het_RA_count </th> <th> hom_RR_count </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 71719 </td> <td align="right"> 2090 </td> <td align="right"> 92239 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 72908 </td> <td align="right"> 1950 </td> <td align="right"> 90665 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 0.59 </td> <td align="right"> 0.41 </td> <td align="right"> 44798 </td> <td align="right"> 64482 </td> <td align="right"> 57929 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 74881 </td> <td align="right"> 2165 </td> <td align="right"> 88837 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 73620 </td> <td align="right"> 2592 </td> <td align="right"> 89160 </td> </tr>
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

Ethnicity Inference
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Sample-Level-QC.md#ethnicity-inference

Genome Similarity
-----------------------------------
https://github.com/deflaux/codelabs/blob/qc-codelab/R/PlatinumGenomes-QC/Sample-Level-QC.md#genome-similarity

Ti/Tv by Depth
-----------------------------------

Heterozygous haplotype
-----------------------------------








