SELECT 
reference_name,
start,
end,
missingness_rate,
FROM (
  SELECT
  reference_name,
  start,
  end,
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
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
    SUM(call.genotype == -1) WITHIN RECORD AS no_calls,
    COUNT(call.genotype) WITHIN RECORD AS all_calls,
    FROM
    [_THE_EXPANDED_TABLE_]
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
missingness_rate > _CUTOFF_
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
