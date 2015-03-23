SELECT
  genotype,
  COUNT(genotype) AS cnt
FROM
  (
  SELECT
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    call.FILTER,
    genotype
  FROM
    (
    SELECT
      reference_name,
      start,
      reference_bases,
      call.FILTER,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotype
    FROM 
      FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs], call.call_set_name)
    WHERE
      call.FILTER = "PASS"
    )
  GROUP EACH BY
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    call.FILTER,
    genotype
  )
GROUP BY
  genotype
