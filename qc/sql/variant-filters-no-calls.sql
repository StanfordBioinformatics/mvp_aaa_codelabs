SELECT
  call.call_set_name AS sample_name,
  call.FILTER AS filter,
  genotype,
  COUNT(call.FILTER) AS count
FROM (
  SELECT
    call.call_set_name,
    call.FILTER,
    alternate_bases,
    GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
  FROM (
    FLATTEN([_THE_TABLE_], call.call_set_name)))
WHERE 
  genotype = '-1,-1'
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  sample_name,
  filter,
  genotype,
ORDER BY
  sample_name,
  count DESC