SELECT
  reference_name,
  VAR_type,
  COUNT(VAR_type) AS Cnt
FROM
  (
  SELECT
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    VAR_type
  FROM
    (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      VAR_type
    FROM
      (
      SELECT
        reference_name,
        start,
        reference_bases,
        alternate_bases,
        IF(LENGTH(reference_bases)=1 AND LENGTH(alternate_bases)=1, "SNV", "INDEL") AS VAR_type,
        call.FILTER
      FROM 
        FLATTEN([gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_vcfs], call.call_set_name)
      WHERE
        call.FILTER = "PASS"
      )
    )
  GROUP EACH BY
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    VAR_type
  )
GROUP BY
  reference_name,
  VAR_type
ORDER BY
  reference_name
