SELECT
sample.id,
(hom_AA_count + het_RA_count) AS all_SNV_count,
hom_AA_count,
het_RA_count,
ROUND((het_RA_count/hom_AA_count)*1000)/1000 AS Het_Hom_ratio
FROM
(
  SELECT
  call.call_set_name AS sample.id,
  SUM(genotype="1/1" AND (LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")))) AS hom_AA_count,
  SUM((genotype="0/1" OR genotype="1/0" OR genotype="1/2") AND (LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")))) AS het_RA_count    
  FROM
  (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotype,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    call.FILTER
    FROM 
    FLATTEN([_THE_EXPANDED_TABLE_], alternate_bases)
    OMIT call IF SOME(call.qc IS NOT NULL)
    HAVING
    cohort_qc IS NULL AND call.FILTER = "PASS"
  ) 
  GROUP BY
  sample.id
)
GROUP BY
sample.id,
hom_AA_count,
het_RA_count,
all_SNV_count,
Het_Hom_ratio
ORDER BY
sample.id;