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
      [_THE_EXPANDED_TABLE_]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
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
