# Count the number of SNPs found in dbSNP with the multi sample variants table.
# A position is considered to be in dbSNP only if there is an rsid in the 'names'
# column.

SELECT
call.call_set_name AS sample_id,
VAR_type,
COUNT(call.call_set_name) AS num_VAR,
COUNT(dbSNP) AS num_VAR_dbSNP
FROM (
  SELECT
  call.call_set_name,
  IF(LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")) , "SNV", "INDEL") AS VAR_type,
  dbSNP
  FROM (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(names) WITHIN record AS dbSNP,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    call.FILTER
    FROM 
    FLATTEN([_THE_EXPANDED_TABLE_], alternate_bases)
    OMIT call IF SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL AND call.FILTER = "PASS"
  )
)
GROUP BY
sample_id,
VAR_type
ORDER BY
sample_id,
VAR_type
