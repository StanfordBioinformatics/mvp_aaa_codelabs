# Comparison of Google Genomics tools vs Academic Tools

The purpose of this codelab is to compare BigQuery queries to academic command line tools for genomic analysis. 

* [Setup](#setup)  
* [Singleton Rate](#singleton-rate)
* [Homozygosity Rate and Inbreeding Coefficient](#homozygosity-rate-and-inbreeding-coefficient)
* [Hardy-Weinberg Equilibrium](#hardy-weingberg-equilibrium)
* [Ti/Tv by Alternate Allele Counts](#titv-by-alternate-allele-counts)
* [Ti/Tv by Depth](#titv-by-depth)
* [Sample-Level Missingness](#sample-level-missingness)
* [Variant-Level Missingness](#variant-level-missingness)
* [Sex Inference](#sex-inference)
* [Heterozygous Haplotype](#heterozygous-haplotype)
* [Ethnicity Inference](#ethnicity-inference)
* [Genome Similarity](#genome-similarity)
* [Genotype Concordance](#genotype-concordance)

## Setup






```r
tableReplacement <- list("_THE_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2",
                          "_THE_EXPANDED_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2",
                         "_REF_TABLE_"="gbsc-gcp-project-mvp:qc_tables.5_genomes_ref_calls_brca1",
                         "_VARIANT_TABLE_"="gbsc-gcp-project-mvp:qc_tables.5_genomes_variants_brca1",
                         "_GENOTYPING_TABLE_"="gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_genotyping_vcfs")
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

constraints <- list("#_WHERE_"="WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499")

googleRepository = "https://raw.githubusercontent.com/googlegenomics/codelabs/master/R/PlatinumGenomes-QC/sql/"
```
 

[Top](#top)

## Singleton Rate
[Singletone Rate in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#singleton-rate)

For each genome, count the number of variants shared by no other member of the cohort.  Too many private calls for a particular individual may indicate a problem.

*For our 5 genome sample set the singleton rate will be very high due to the small number of genomes.*


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
          [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
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
<!-- Mon May  4 15:50:24 2015 -->
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

Command line execution
```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --singletons --out singletons
```


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

How many singletons were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult)
```

```
## Joining by: c("CHROM", "POS", "SINGLETON_DOUBLETON", "INDV")
```

```r
nrow(onlyBQ)
```

```
## [1] 0
```


```
## Error in matrix(unlist(value, recursive = FALSE, use.names = FALSE), nrow = nr, : length of 'dimnames' [2] not equal to array extent
```

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
<!-- Mon May  4 15:50:24 2015 -->
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
  <tr> <td> chr17 </td> <td align="right"> 41230105.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41229760.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226740.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226739.00 </td> <td> D </td> <td> T </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226738.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226737.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41226736.00 </td> <td> D </td> <td> G </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41214210.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41242075.00 </td> <td> D </td> <td> A </td> <td> LP6005038-DNA_B02 </td> </tr>
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
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
WHERE
  reference_name = 'chr17'
HAVING
  start = 41271293 OR start = 41267518 OR start = 41256102 OR start = 41256101 OR start = 41256100 OR start = 41256086 OR start = 41256085 OR start = 41256084 OR start = 41256083 OR start = 41256082 OR start = 41256081 OR start = 41256080 OR start = 41256079 OR start = 41256078 OR start = 41256077 OR start = 41256076 OR start = 41256075 OR start = 41252696 OR start = 41252647 OR start = 41252634 OR start = 41249363 OR start = 41242077 OR start = 41242076 OR start = 41230105 OR start = 41229760 OR start = 41226740 OR start = 41226739 OR start = 41226738 OR start = 41226737 OR start = 41226736 OR start = 41214210 OR start = 41242075 OR start = 41214209 OR start = 41202634
ORDER BY
  start,
  end,
  call.call_set_name
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:27 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> genotype </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th>  </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41252696 </td> <td align="right"> 41252755 </td> <td> T </td> <td>  </td> <td align="right"> 31.21 </td> <td>  </td> </tr>
   </table>

It appears that they correspond either to:

* A reference-matching block, so not actually a singleton and just perhaps violating an assumption in the vcftools code.  Most of the start sites queried return nothing because they are withing reference matching blocks.
* Or a non-singleon variant, perhaps due to a problem in converting the gVCF data to all-positions VCF via gvcftools?

[Top](#top)

## Homozygosity Rate and Inbreeding Coefficient
[Homozygosity Rate and Inbreeding Coefficient in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#homozygosity-rate-and-inbreeding-coefficient)

For each genome, compare the expected and observed rates of homozygosity.


```r
query <- paste(googleRepository, "homozygous-variants.sql", sep="") 
result <- DisplayAndDispatchQuery(query,
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
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
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
<!-- Mon May  4 15:50:30 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> O_HOM </th> <th> E_HOM </th> <th> N_SITES </th> <th> F </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right">  16 </td> <td align="right"> 79.87 </td> <td align="right"> 136 </td> <td align="right"> -1.14 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right">  14 </td> <td align="right"> 68.87 </td> <td align="right"> 135 </td> <td align="right"> -0.83 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right">  15 </td> <td align="right"> 68.23 </td> <td align="right"> 134 </td> <td align="right"> -0.81 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right">  18 </td> <td align="right"> 80.00 </td> <td align="right"> 136 </td> <td align="right"> -1.11 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 137 </td> <td align="right"> 81.97 </td> <td align="right"> 140 </td> <td align="right"> 0.95 </td> </tr>
   </table>

Let's compare to what was found using vcftools.

Command line execution.

```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --het --out het_calls
```


```r
expectedResult <- read.table("./data/brca1/het_calls.het", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:30 2015 -->
<table border=1>
<tr> <th> E_HOM.x </th> <th> E_HOM.y </th> <th> F.x </th> <th> F.y </th> <th> INDV </th> <th> N_SITES.x </th> <th> N_SITES.y </th> <th> O_HOM.x </th> <th> O_HOM.y </th>  </tr>
  <tr> <td align="right"> 65.50 </td> <td align="right"> 79.87 </td> <td align="right"> -0.85 </td> <td align="right"> -1.14 </td> <td> LP6005038-DNA_A01 </td> <td align="right"> 143 </td> <td align="right"> 136 </td> <td align="right">   0 </td> <td align="right">  16 </td> </tr>
  <tr> <td align="right"> 66.40 </td> <td align="right"> 68.87 </td> <td align="right"> -0.70 </td> <td align="right"> -0.83 </td> <td> LP6005038-DNA_A02 </td> <td align="right"> 150 </td> <td align="right"> 135 </td> <td align="right">   8 </td> <td align="right">  14 </td> </tr>
  <tr> <td align="right"> 67.60 </td> <td align="right"> 68.23 </td> <td align="right"> -0.72 </td> <td align="right"> -0.81 </td> <td> LP6005038-DNA_A03 </td> <td align="right"> 149 </td> <td align="right"> 134 </td> <td align="right">   9 </td> <td align="right">  15 </td> </tr>
  <tr> <td align="right"> 64.00 </td> <td align="right"> 80.00 </td> <td align="right"> -0.87 </td> <td align="right"> -1.11 </td> <td> LP6005038-DNA_B01 </td> <td align="right"> 135 </td> <td align="right"> 136 </td> <td align="right">   2 </td> <td align="right">  18 </td> </tr>
  <tr> <td align="right"> 61.80 </td> <td align="right"> 81.97 </td> <td align="right"> 0.92 </td> <td align="right"> 0.95 </td> <td> LP6005038-DNA_B02 </td> <td align="right"> 133 </td> <td align="right"> 140 </td> <td align="right"> 127 </td> <td align="right"> 137 </td> </tr>
   </table>

The logic in the query looks similar to vcftools [output_het method](http://sourceforge.net/p/vcftools/code/HEAD/tree/trunk/cpp/variant_file_output.cpp#l165) but there is clearly a difference.  

TODO: investigate the difference further.

[Top](#top)

Hardy-Weinberg Equilibrium
-----------------------------------
[Hardy-Weinberg Equiibrium in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Variant-Level-QC.md#hardy-weinberg-equilibrium)


```r
query <- paste(googleRepository, "hardy-weinberg.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
# The following query computes the Hardy-Weinberg equilibrium for variants.
SELECT
  reference_name,
  start,
  reference_bases,
  alternate_bases,
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
      reference_name,
      start,
      reference_bases,
      alternate_bases,
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
      reference_name,
      start,
      reference_bases,
      alternate_bases,
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
        [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
      # Optionally add a clause here to limit the query to a particular
      # region of the genome.
      WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        )))
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
```
Number of rows returned by this query: 185.

Displaying the first few results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:33 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alternate_bases </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41277186 </td> <td> G </td> <td> C </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> <td align="right"> 2.22 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230104 </td> <td> CT </td> <td> C </td> <td align="right">   0 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right"> 4.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230227 </td> <td> G </td> <td> A </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> <td align="right"> 2.22 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230335 </td> <td> A </td> <td> G </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> <td align="right"> 2.22 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230375 </td> <td> A </td> <td> G </td> <td align="right">   1 </td> <td align="right">   4 </td> <td align="right">   0 </td> <td align="right"> 1.80 </td> <td align="right"> 2.40 </td> <td align="right"> 0.80 </td> <td align="right"> 2.22 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41230523 </td> <td> T </td> <td> G </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> <td align="right"> 0.31 </td> <td> FALSE </td> </tr>
   </table>

Compare to [brca1.hwe](./data/brca1/hardy.hwe) (see the [vcftools command line](./data/hwe/brca1.log) used to create this file).

Command line
```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --hardy --out hardy
```


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
names(expectedResult)[1] = 'reference_name'
names(expectedResult)[2] = 'start'
```

How many results do the two results have in common?

```r
nrow(inner_join(result, expectedResult, by=c("reference_name", "start", "OBS_HOM1", "OBS_HET", "OBS_HOM2")))
```

```
## [1] 160
```

How many result were only identified by BigQuery?

```r
onlyBQ <- anti_join(result, expectedResult, , by=c("reference_name", "start", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
nrow(onlyBQ)
```

```
## [1] 25
```

Let's take a look at the first few.
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:33 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> reference_bases </th> <th> alternate_bases </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th> <th> ChiSq </th> <th> PVALUE_SIG </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944 </td> <td> T </td> <td> C </td> <td align="right">   4 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 4.05 </td> <td align="right"> 0.90 </td> <td align="right"> 0.05 </td> <td align="right"> 0.06 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201363 </td> <td> T </td> <td> C </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> <td align="right"> 0.31 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> C </td> <td> CACA </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> <td align="right"> 1.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41211485 </td> <td> C </td> <td> CACAACA </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> <td align="right"> 1.00 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213759 </td> <td> C </td> <td> T </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> <td align="right"> 0.31 </td> <td> FALSE </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216205 </td> <td> T </td> <td> C </td> <td align="right">   3 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 3.20 </td> <td align="right"> 1.60 </td> <td align="right"> 0.20 </td> <td align="right"> 0.31 </td> <td> FALSE </td> </tr>
   </table>

How many results were only identified by vcftools?

```r
onlyVcftools <- anti_join(expectedResult, result, , by=c("reference_name", "start", "OBS_HOM1", "OBS_HET", "OBS_HOM2"))
nrow(onlyVcftools)
```

```
## [1] 16
```

Let's take a look at the first few.
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:33 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> ChiSq </th> <th> P </th> <th> OBS_HOM1 </th> <th> OBS_HET </th> <th> OBS_HOM2 </th> <th> E_HOM1 </th> <th> E_HET </th> <th> E_HOM2 </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196944.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41201363.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41213759.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41216205.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41220771.00 </td> <td align="right"> 1.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right">   0 </td> <td align="right"> 0.25 </td> <td align="right"> 0.50 </td> <td align="right"> 0.25 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41222098.00 </td> <td align="right"> 2.00 </td> <td align="right"> 1.00 </td> <td align="right">   0 </td> <td align="right">   2 </td> <td align="right">   0 </td> <td align="right"> 0.50 </td> <td align="right"> 1.00 </td> <td align="right"> 0.50 </td> </tr>
   </table>

Retrieving the gVCF data for the results identified only by vcftools:

```r
having <- list("_HAVING_" = paste("start <= ", onlyVcftools$start,
                "AND",
                "end >= ", onlyVcftools$start+1,
                collapse=" OR "))
query <- paste(googleRepository, "examine-data.sql", sep="")
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement,
                                                 having))
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
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
WHERE
  reference_name = 'chr17'
HAVING
  start <=  41258790 AND end >=  41258791 OR start <=  41252694 AND end >=  41252695 OR start <=  41249362 AND end >=  41249363 OR start <=  41248903 AND end >=  41248904 OR start <=  41247623 AND end >=  41247624 OR start <=  41247054 AND end >=  41247055 OR start <=  41245470 AND end >=  41245471 OR start <=  41241567 AND end >=  41241568 OR start <=  41239979 AND end >=  41239980 OR start <=  41251930 AND end >=  41251931 OR start <=  41222098 AND end >=  41222099 OR start <=  41220771 AND end >=  41220772 OR start <=  41216205 AND end >=  41216206 OR start <=  41213759 AND end >=  41213760 OR start <=  41201363 AND end >=  41201364 OR start <=  41196944 AND end >=  41196945
ORDER BY
  start,
  end,
  call.call_set_name
```

The first few:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:36 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> genotype </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> reference_bases </th> <th> alternate_bases </th> <th> quality </th> <th> filter </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196408 </td> <td align="right"> 41197273 </td> <td> G </td> <td>  </td> <td align="right"> 61.23 </td> <td>  </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196408 </td> <td align="right"> 41197273 </td> <td> G </td> <td>  </td> <td align="right"> 61.23 </td> <td>  </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196408 </td> <td align="right"> 41197273 </td> <td> G </td> <td>  </td> <td align="right"> 61.23 </td> <td>  </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41196408 </td> <td align="right"> 41197273 </td> <td> G </td> <td>  </td> <td align="right"> 61.23 </td> <td>  </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td> 0,1 </td> <td> chr17 </td> <td align="right"> 41196944 </td> <td align="right"> 41196945 </td> <td> T </td> <td> C </td> <td align="right"> 315.77 </td> <td> PASS </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> 0,0 </td> <td> chr17 </td> <td align="right"> 41200109 </td> <td align="right"> 41201701 </td> <td> A </td> <td>  </td> <td align="right"> 40.23 </td> <td>  </td> </tr>
   </table>

[Top](#top)

Ti/Tv by Alternate Allele Counts
-----------------------------------
[Ti/Tv by Alternate Allele Counts in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Variant-Level-QC.md#titv-by-alternate-allele-counts)


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
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
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
<!-- Mon May  4 15:50:38 2015 -->
<table border=1>
<tr> <th> transitions </th> <th> transversions </th> <th> titv </th> <th> alternate_allele_count </th>  </tr>
  <tr> <td align="right">   0 </td> <td align="right">   1 </td> <td align="right"> 0.00 </td> <td align="right">  10 </td> </tr>
  <tr> <td align="right">  81 </td> <td align="right">  29 </td> <td align="right"> 2.79 </td> <td align="right">   4 </td> </tr>
  <tr> <td align="right">  11 </td> <td align="right">   3 </td> <td align="right"> 3.67 </td> <td align="right">   2 </td> </tr>
  <tr> <td align="right">   8 </td> <td align="right">   5 </td> <td align="right"> 1.60 </td> <td align="right">   1 </td> </tr>
  <tr> <td align="right">   3 </td> <td align="right">   0 </td> <td align="right">  </td> <td align="right">   0 </td> </tr>
   </table>
    
Compare to [output](./data/brca1/TsTv-by-count.TsTv.count) from vcftools.

Command line
```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --TsTv-by-count --out TsTv-by-coun
```


```r
expectedResult <- read.table("./data/brca1/TsTv-by-count.TsTv.count", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:38 2015 -->
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

The results for BigQuery and vcftools are identical except for the case where there are 10 alternate alleles.  With only 5 genomes, 10 sites with an alternate is the maximum.  vcftools may ignore sites where all calls are the alternate allele.

[Top](#top)

Ti/Tv by Depth
-----------------------------------

```r
query <- "./sql/ti-tv-by-depth.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement))
```

```
SELECT
  call.call_set_name,
  (transitions/transversions) AS titv_ratio,
  average_depth,
FROM (
  SELECT
    call.call_set_name,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    ROUND(AVG(call.DP)) AS average_depth,
  FROM (

    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.DP
    FROM (
      SELECT
        call.call_set_name,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        call.genotype,
        call.DP,
      FROM
        [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
      # Optionally add clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_  
      )
    WHERE
      call.DP is not null
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
    GROUP BY 
      call.call_set_name,
      call.DP,)
WHERE
  transversions > 0
GROUP BY
  call.call_set_name,
  titv_ratio,
  average_depth,
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:41 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> titv_ratio </th> <th> average_depth </th>  </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 2.09 </td> <td align="right"> 39.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 2.20 </td> <td align="right"> 28.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 2.20 </td> <td align="right"> 30.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 2.12 </td> <td align="right"> 30.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 2.15 </td> <td align="right"> 35.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 2.16 </td> <td align="right"> 33.00 </td> </tr>
   </table>


```r
ggplot(result, aes(x=average_depth, y=titv_ratio, color=call_call_set_name)) + 
  geom_point() +
  #geom_smooth(se=FALSE,    # Don't add shaded confidence region
  #              fullrange=T) +
  ggtitle("Ti/Tv Ratio By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv")
```

<img src="figure/titv-by-depth-1.png" title="plot of chunk titv-by-depth" alt="plot of chunk titv-by-depth" style="display: block; margin: auto;" />


TODO: Find academic tool to compare against.

[Top](#top)

Sample-Level Missingness
-----------------------------------
[Sample-Level Missingness in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#missingness-rate)

See [this document](./Sample_Missingness_Preparation.md) for details on how tables were created for this section.


```r
query <- "./sql/sample-missingness-missing-calls.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
SELECT 
  bc.sample_name AS sample_name,
  nc.no_call_count as no_call_count,
  nc.no_call_count/bc.base_count AS missingness,
FROM (
  SELECT 
    sample_name,
    COUNT(no_call) as no_call_count,
  FROM(
    SELECT
      ref.sample_name AS sample_name,
      ref.chromosome AS chromosome,
      ref.start AS start,
      ref.end AS end,
      ref.reference_bases AS reference_bases,
      True AS no_call
    FROM
      (SELECT 
        *
       FROM
        [gbsc-gcp-project-mvp:qc_tables.5_genomes_ref_calls_brca1] ) AS ref
    LEFT OUTER JOIN (
      SELECT
        *
      FROM
        [gbsc-gcp-project-mvp:qc_tables.5_genomes_variants_brca1] ) AS var
    ON
      ref.sample_name = var.sample_name
      AND ref.chromosome = var.chromosome
      AND ref.start = var.start
    WHERE
      ref.is_ref is False
      AND var.is_variant_call is NULL )
  GROUP BY
    sample_name ) AS nc
JOIN (
  SELECT
    sample_name,
    COUNT(start) AS base_count
  FROM 
    [gbsc-gcp-project-mvp:qc_tables.5_genomes_ref_calls_brca1]
  GROUP BY
  sample_name ) as bc
ON
  bc.sample_name = nc.sample_name;
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:44 2015 -->
<table border=1>
<tr> <th> sample_name </th> <th> no_call_count </th> <th> missingness </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 286 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 115 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 332 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 228 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 639 </td> <td align="right"> 0.01 </td> </tr>
   </table>
    
Compare to [output](./data/brca1/missingness.imiss) from vcftools.

Command line
```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --missing-indv --out missingness
```


```r
expectedResult <- read.table("./data/brca1/missingness.imiss", header=TRUE)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:44 2015 -->
<table border=1>
<tr> <th> INDV </th> <th> N_DATA </th> <th> N_GENOTYPES_FILTERED </th> <th> N_MISS </th> <th> F_MISS </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 274 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 107 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 311 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right"> 205 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 81239 </td> <td align="right">   0 </td> <td align="right">  73 </td> <td align="right"> 0.00 </td> </tr>
   </table>

These results from BigQuery and vcftools are expected to be different.  Sample-level missingness in BigQuery is calculated by joining all positions in the genome against the reference genome to identify all sites that were not called by sequencing.  Vcftools only looks at sites that were called by at least one sample in the group.

[Top](#top)

Variant-Level Missingness
-----------------------------------
[Variant-Level Missingness in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Variant-Level-QC.md#missingness-rate)

```r
query <- "./sql/variant-level-missingness-missing-calls.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints))
```

```
#SELECT  FROM [va_aaa_pilot_data.5_genome_test_vcfs_no_calls] LIMIT 1000
# Compute the ratio no-calls for each variant.
SELECT
  reference_name,
  start,
  END,
  reference_bases,
  alternate_bases,
  no_calls,
  all_calls,
  (sample_count*2 - all_calls) AS missing_calls,
  ((no_calls + (sample_count*2 - all_calls))/(all_calls + (sample_count*2 - all_calls))) AS missingness_rate,
  
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
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    WHERE 
                    reference_name = 'chr17'
                    AND start BETWEEN 41196311
                    AND 41277499
  ) as calls
CROSS JOIN (
  SELECT 
    COUNT(sample_name) AS sample_count
  FROM(
    SELECT
      call.call_set_name AS sample_name, 
    FROM [va_aaa_pilot_data.5_genome_test_vcfs_no_calls] 
    GROUP BY sample_name)) as count
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
```

```
Error: Not found: Table gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls

 notFound. Not found: Table gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_no_calls
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:45 2015 -->
<table border=1>
<tr> <th> sample_name </th> <th> no_call_count </th> <th> missingness </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 286 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 115 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 332 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 228 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 639 </td> <td align="right"> 0.01 </td> </tr>
   </table>
    
    
Compare to [output](./data/brca1/missing_sites.lmiss) from vcftools.

Command line
```
vcftools --gzvcf brca1.merged.sample_set.vcf.gz --missing_site --out missing_sites
```


```r
expectedResult <- read.table("./data/brca1/missing_sites.lmiss", header=TRUE)
```



```r
expectedResult <- mutate(expectedResult, POS = POS - 1)
n = c("reference_name", "start", "expected_count", "filtered", "expected_missing_count", "expected_rate")
names(expectedResult) = n
```


```r
print(xtable(head(expectedResult)), type="html", include.rownames=F)
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:50:46 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> expected_count </th> <th> filtered </th> <th> expected_missing_count </th> <th> expected_rate </th>  </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196311.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196312.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196313.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196314.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196315.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
  <tr> <td> chr17 </td> <td align="right"> 41196316.00 </td> <td align="right">  10 </td> <td align="right">   0 </td> <td align="right">   0 </td> <td align="right"> 0.00 </td> </tr>
   </table>

Let's join them and look at the differences

```r
joinedResult = data.table(inner_join(result, expectedResult, by=c("reference_name", "start")))
```

```
## Error in eval(expr, envir, enclos): could not find function "data.table"
```

```r
matches = joinedResult[missingness_rate == expected_rate]
```

```
## Error in `[.data.frame`(joinedResult, missingness_rate == expected_rate): object 'missingness_rate' not found
```

```r
differences = joinedResult[missingness_rate != expected_rate]
```

```
## Error in `[.data.frame`(joinedResult, missingness_rate != expected_rate): object 'missingness_rate' not found
```

How many do they have in common?

```r
nrow(matches)
```

```
## Error in nrow(matches): object 'matches' not found
```

How many are different?

```r
nrow(differences)
```

```
## Error in nrow(differences): object 'differences' not found
```

Here we can see which rows where different

```
## Error in xtable(differences): object 'differences' not found
```

The differences appear to be mostly form positions with insertions or deletions.

TODO: Look at the differences more closely.

[Top](#top)

Sex Inference
-----------------------------------
[Sex Inference in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#sex-inference)


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
    [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
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
<!-- Mon May  4 15:50:49 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> perct_het_alt_in_snvs </th> <th> perct_hom_alt_in_snvs </th> <th> hom_AA_count </th> <th> het_RA_count </th> <th> hom_RR_count </th>  </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 71719 </td> <td align="right"> 2090 </td> <td align="right"> 91043 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 72908 </td> <td align="right"> 1950 </td> <td align="right"> 89501 </td> </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 0.59 </td> <td align="right"> 0.41 </td> <td align="right"> 44798 </td> <td align="right"> 64482 </td> <td align="right"> 56682 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 74881 </td> <td align="right"> 2165 </td> <td align="right"> 87670 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 0.03 </td> <td align="right"> 0.97 </td> <td align="right"> 73620 </td> <td align="right"> 2592 </td> <td align="right"> 88003 </td> </tr>
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

TODO: Compare to Plink

Command line
```
plink --vcf ~/scratch/chrX.vcf.gz --make-bed --split-x 'hg19' --out ~/scratch/chrX-split
plink --bfile ~/scratch/chrX-split --check-sex
```

[Top](#top)

Heterozygous Haplotype
-----------------------------------
[Heterozygous Haplotype in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Variant-Level-QC.md#heterozygous-haplotype)

```r
query <- paste(googleRepository, "sex-chromosome-heterozygous-haplotypes.sql", sep="")
male_sample_ids = paste("'", paste(sampleInfo[sampleInfo$gender == 'Male',]$call_call_set_name, collapse=",'", "'", sep=''))
male_sample_ids = gsub(" ", "", male_sample_ids, fixed = TRUE)
male_sub = list("_MALE_SAMPLE_IDS_" = male_sample_ids)
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, constraints, male_sub))
```

```
# Retrieve heterozygous haplotype calls on chromosomes X and Y.
SELECT
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
FROM
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
WHERE
  reference_name IN ('chrX', 'chrY')
OMIT
  call if (2 > COUNT(call.genotype))
  OR EVERY(call.genotype <= 0)
  OR EVERY(call.genotype = 1)
HAVING call.call_set_name IN ('')
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
```


```
Error in UseMethod("xtable"): no applicable method for 'xtable' applied to an object of class "NULL"
```

TODO Compare to command line

[Top](#top)

Ethnicity Inference
-----------------------------------
[Ethnicity Inference in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#ethnicity-inference)

Status: Waiting for PCA release from Elmer.  See [here](https://github.com/elmer-garduno/spark-examples/tree/multiple_dataset_pca) for status.

[Top](#top)

Genome Similarity
-----------------------------------
[Genome Similarity in Google Genomics](https://github.com/googlegenomics/codelabs/blob/master/R/PlatinumGenomes-QC/Sample-Level-QC.md#genome-similarity)

Status: Installed apache maven and built Google Genomics dataflow package.  Ran into errors when trying to run dataflow job.

```
java -cp target/google-genomics-dataflow-v1beta2-0.2-SNAPSHOT.jar   com/google/cloud/genomics/dataflow/pipelines/IdentityByState   --project=my-project-id   --output=gs://my-bucket/localtest.txt   --genomicsSecretsFile=client_secrets.json
Error: Could not find or load main class com.google.cloud.genomics.dataflow.pipelines.IdentityByState
```

Command line
```
plink --genome --vcf chr22-merged.vcf.gz
```

[Top](#top)

Genotype Concordance
-----------------------------------

Calculate the concordance between sequencing data and genotyping data.


```r
query <- "./sql/genotyping-concordance.sql"
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement))
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
          [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
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
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_genotyping_vcfs]
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
Running query:   RUNNING  2.6s
Running query:   RUNNING  3.2s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.7s
Running query:   RUNNING  6.3s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.5s
Running query:   RUNNING  8.1s
Running query:   RUNNING  8.8s
Running query:   RUNNING  9.4s
Running query:   RUNNING 10.0s
Running query:   RUNNING 10.6s
Running query:   RUNNING 11.2s
Running query:   RUNNING 11.9s
Running query:   RUNNING 12.5s
Running query:   RUNNING 13.1s
Running query:   RUNNING 13.7s
Running query:   RUNNING 14.3s
Running query:   RUNNING 14.9s
Running query:   RUNNING 15.6s
Running query:   RUNNING 16.2s
Running query:   RUNNING 16.8s
Running query:   RUNNING 17.4s
Running query:   RUNNING 18.0s
Running query:   RUNNING 18.7s
Running query:   RUNNING 19.3s
Running query:   RUNNING 19.9s
Running query:   RUNNING 20.5s
Running query:   RUNNING 21.1s
Running query:   RUNNING 21.8s
Running query:   RUNNING 22.4s
Running query:   RUNNING 23.0s
Running query:   RUNNING 23.6s
Running query:   RUNNING 24.2s
Running query:   RUNNING 24.8s
Running query:   RUNNING 25.4s
Running query:   RUNNING 26.1s
Running query:   RUNNING 26.7s
Running query:   RUNNING 27.3s
Running query:   RUNNING 27.9s
Running query:   RUNNING 28.5s
Running query:   RUNNING 29.2s
Running query:   RUNNING 29.8s
Running query:   RUNNING 30.4s
Running query:   RUNNING 31.0s
Running query:   RUNNING 31.6s
Running query:   RUNNING 32.3s
Running query:   RUNNING 32.9s
Running query:   RUNNING 33.5s
Running query:   RUNNING 34.1s
Running query:   RUNNING 34.8s
Running query:   RUNNING 35.4s
Running query:   RUNNING 36.0s
Running query:   RUNNING 36.6s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.8s
Running query:   RUNNING 38.5s
Running query:   RUNNING 39.1s
Running query:   RUNNING 39.7s
Running query:   RUNNING 40.3s
Running query:   RUNNING 40.9s
Running query:   RUNNING 41.6s
Running query:   RUNNING 42.2s
Running query:   RUNNING 42.8s
Running query:   RUNNING 43.4s
```

<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Mon May  4 15:51:37 2015 -->
<table border=1>
<tr> <th> sample_id </th> <th> calls_in_common </th> <th> identical_calls </th> <th> concordance </th>  </tr>
  <tr> <td> LP6005038-DNA_A03 </td> <td align="right"> 2199423 </td> <td align="right"> 2177458 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td align="right"> 2200890 </td> <td align="right"> 2178871 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td align="right"> 2200420 </td> <td align="right"> 2178471 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td align="right"> 2200848 </td> <td align="right"> 2178903 </td> <td align="right"> 0.99 </td> </tr>
  <tr> <td> LP6005038-DNA_B02 </td> <td align="right"> 2201742 </td> <td align="right"> 2179131 </td> <td align="right"> 0.99 </td> </tr>
   </table>






