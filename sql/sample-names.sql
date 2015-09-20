# Select only the sample names from a table

SELECT
  call.call_set_name AS sample_name,
FROM 
  [va_aaa_pilot_data.5_genome_test_gvcfs_2]
GROUP BY
  sample_name,
ORDER BY
  sample_name,
