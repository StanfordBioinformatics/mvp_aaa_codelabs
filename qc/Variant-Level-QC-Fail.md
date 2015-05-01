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

# Variant-Level QC





* [Missingness Rate](#missingness-rate)
* [Ti/Tv by Depth](#titv-by-depth)



```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs")
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

# To run this against other public data, source in one of the dataset helpers.  For example:
# source("./rHelpers/pgpCGIOnlyDataset.R")
```

## Missingness Rate

For each variant, compute the missingness rate.  This query can be used to identify variants with a poor call rate.


```r
cutoff = list("_CUTOFF_"="0.9")
result <- DisplayAndDispatchQuery("./sql/variant-level-missingness-fail.sql",
                                  project=project,
                                  replacements=c(cutoff,
                                                 queryReplacements))
```

```
# Compute the ratio no-calls for each variant.
SELECT 
reference_name,
start,
END,
missingness_rate,
FROM (
  SELECT
  reference_name,
  start,
  END,
  reference_bases,
  alternate_bases,
  no_calls,
  all_calls,
  (no_calls/all_calls) AS no_call_rate,
  1 - (all_calls-no_calls)/sample_count AS missingness_rate,
  sample_count
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
    [va_aaa_pilot_data.all_genomes_expanded_vcfs]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
  ) as.g
  CROSS JOIN (
    SELECT
    COUNT(call.call_set_name) AS sample_count
    FROM (
      SELECT 
      call.call_set_name
      FROM
      [va_aaa_pilot_data.all_genomes_gvcfs]
      GROUP BY 
      call.call_set_name)) AS count )
WHERE
missingness_rate > 0.9
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_

Retrieving data:  3.3s
Retrieving data:  4.5s
Retrieving data:  5.7s
Retrieving data:  7.3s
Retrieving data:  8.7s
Retrieving data: 10.0s
Retrieving data: 11.4s
Retrieving data: 12.6s
```
Number of rows returned by this query: **100000**.

Displaying the first few rows of the dataframe of results:
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May  1 02:48:22 2015 -->
<table border=1>
<tr> <th> reference_name </th> <th> start </th> <th> END </th> <th> missingness_rate </th>  </tr>
  <tr> <td> chr3 </td> <td align="right"> 73056435 </td> <td align="right"> 73056437 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr3 </td> <td align="right"> 73908346 </td> <td align="right"> 73908347 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr3 </td> <td align="right"> 73908354 </td> <td align="right"> 73908355 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr3 </td> <td align="right"> 73908357 </td> <td align="right"> 73908358 </td> <td align="right"> 0.97 </td> </tr>
  <tr> <td> chr3 </td> <td align="right"> 73908357 </td> <td align="right"> 73908358 </td> <td align="right"> 1.00 </td> </tr>
  <tr> <td> chr3 </td> <td align="right"> 73908357 </td> <td align="right"> 73908358 </td> <td align="right"> 1.00 </td> </tr>
   </table>


## Ti/Tv by Depth


```r
query <- "./sql/variant-depth-fail.sql"
max = 70
min = 8
cutoffs = list("_MAX_VALUE_" = max,
               "_MIN_VALUE_" = min)
result <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement, cutoffs))
```

```
SELECT
  call.call_set_name,
  reference_name,
  start,
  end,
  GROUP_CONCAT(STRING(call.DP)) WITHIN call AS call.DP,
FROM
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
WHERE
  call.DP > 70
  OR call.DP < 8
Retrieving data:  3.1s
Retrieving data:  4.8s
Retrieving data:  6.4s
Retrieving data:  7.7s
Retrieving data:  9.1s
Retrieving data: 10.2s
Retrieving data: 11.8s
Retrieving data: 13.2s
```

Number of rows returned by this query: **100000**.

First few results
<!-- html table generated in R 3.1.2 by xtable 1.7-4 package -->
<!-- Fri May  1 02:48:38 2015 -->
<table border=1>
<tr> <th> call_call_set_name </th> <th> reference_name </th> <th> start </th> <th> end </th> <th> call_DP </th>  </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td> chr10 </td> <td align="right"> 1102709 </td> <td align="right"> 1102710 </td> <td> 7 </td> </tr>
  <tr> <td> LP6005038-DNA_A02 </td> <td> chr10 </td> <td align="right"> 1102795 </td> <td align="right"> 1102796 </td> <td> 7 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td> chr10 </td> <td align="right"> 117966666 </td> <td align="right"> 117966667 </td> <td> 3 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td> chr10 </td> <td align="right"> 117966668 </td> <td align="right"> 117966669 </td> <td> 3 </td> </tr>
  <tr> <td> LP6005038-DNA_A01 </td> <td> chr10 </td> <td align="right"> 117966669 </td> <td align="right"> 117966670 </td> <td> 3 </td> </tr>
  <tr> <td> LP6005038-DNA_B01 </td> <td> chr10 </td> <td align="right"> 12463479 </td> <td align="right"> 12463480 </td> <td> 85 </td> </tr>
   </table>


