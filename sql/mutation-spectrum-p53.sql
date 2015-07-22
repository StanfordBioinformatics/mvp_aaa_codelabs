SELECT
call.call_set_name,
reference_name,
window,
window * 1000 AS window_start,
((window * 1000) + 999) AS window_end,
MIN(start) AS min_variant_position,
MAX(start) AS max_variant_position,
COUNT(call.call_set_name) AS num_variants_in_window
FROM (
  SELECT
  call.call_set_name,
  reference_name,
  start,
  end,
  INTEGER(FLOOR(start / 1000)) AS window,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
  GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
  GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
  call.FILTER
  FROM
  [va_aaa_pilot_data.multi_sample_variants_full_qc]
  OMIT 
  call IF EVERY(call.FILTER != "PASS")
  OR SOME(call.QC IS NOT NULL)
  HAVING
  cohort_qc IS NULL
  AND reference_name = "chr17"
  AND ((start <= 7661779 AND end >= 7687550) OR (start >= 7661779 AND start <= 7687550))
)
GROUP BY
call.call_set_name,
reference_name,
window,
window_start,
window_end,
ORDER BY
window ASC,
call.call_set_name ASC;