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
  alternates,
  VAR_type
  FROM
  (
    SELECT
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(QC) WITHIN RECORD AS qc,
    GROUP_CONCAT(call.QC) WITHIN CALL AS call_qc,
    IF(LENGTH(reference_bases)=1 AND LENGTH(alternate_bases)=1, "SNV", "INDEL") AS VAR_type
    FROM 
    [_THE_EXPANDED_TABLE_]
    OMIT call IF EVERY(call.genotype <= 0) 
      OR SOME(call.qc IS NOT NULL)
    HAVING QC IS NULL
  )
  GROUP EACH BY
  reference_name,
  start,
  reference_bases,
  alternates,
  VAR_type
)
GROUP BY
reference_name,
VAR_type
ORDER BY
reference_name,
VAR_type
