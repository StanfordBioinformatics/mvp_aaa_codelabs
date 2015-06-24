SELECT
  call.call_set_name AS sample_name,
  genotype,
  COUNT(genotype) AS count
FROM (
  SELECT
    call.call_set_name,
    alternate_bases,
    GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
  FROM (
    FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_3], call.call_set_name))
)
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  sample_name,
  genotype,
ORDER BY
  sample_name,
  count DESC