SELECT
call.call_set_name AS sample_id,
COUNT(call.call_set_name) AS private_SNVs_count
FROM
(
  SELECT
  call.call_set_name,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples_with_VARs,
  GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
  GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
  FROM
  [_THE_EXPANDED_TABLE_]
  OMIT 
  call IF EVERY(call.FILTER != "PASS") 
  OR SOME(call.QC IS NOT NULL)
  HAVING
  num_samples_with_VARs = 1
  AND cohort_qc IS NULL
)
GROUP EACH BY
sample_id
ORDER BY
sample_id ASC;