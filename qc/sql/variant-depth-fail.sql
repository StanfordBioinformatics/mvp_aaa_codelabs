SELECT
  call.call_set_name as sample_id,
  reference_name,
  start,
  end,
  GROUP_CONCAT(STRING(call.DP)) WITHIN call AS call.DP,
FROM
  [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs_2]
WHERE
  call.DP > _MAX_VALUE_
  OR call.DP < _MIN_VALUE_
