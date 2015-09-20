# Calculate quality metrics (MQ, MQ0, QUAL) for each sample.

SELECT
  call.call_set_name AS sample_name,
  MIN(call.QUAL) AS min_QUAL,
  MIN(MQ) AS min_MQ,
  MAX(MQ0) AS max_MQ0,
  COUNT(call.call_set_name) AS count
FROM (
  SELECT
    call.call_set_name,
    alternate_bases,
    call.QUAL,
    MQ,
    MQ0,
    GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
  FROM (
    FLATTEN([_THE_TABLE_], call.call_set_name)))
WHERE 
  genotype != '-1,-1'
OMIT RECORD IF EVERY(alternate_bases IS NOT NULL)
GROUP BY
  sample_name,
ORDER BY
  sample_name,
  count DESC
