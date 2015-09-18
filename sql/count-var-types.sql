# Count the number of variant types (SNVs or INDELs) for 5 genome test genome calls table.

SELECT
  "5 test genomes" AS sample_id,
  VAR_type,
  COUNT(VAR_type) AS Cnt
FROM (
  SELECT
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    VAR_type
  FROM (
    SELECT
      reference_name,
      start,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      IF(LENGTH(reference_bases)=1 AND LENGTH(alternate_bases)=1, "SNV", "INDEL") AS VAR_type,
    FROM (
      SELECT
        reference_name,
        start,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        call.FILTER,
      FROM 
        [gbsc-gcp-project-mvp:va_aaa_pilot_data.5_genome_test_gvcfs_2]
    )
    WHERE
      call.FILTER = "PASS"
      #_AND_
  )
  GROUP EACH BY
    reference_name,
    start,
    reference_bases,
    alternate_bases,
    VAR_type
)
GROUP BY
  sample_id,
  VAR_type
