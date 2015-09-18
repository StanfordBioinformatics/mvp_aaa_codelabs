# Select variants within ACMG genes.  This query joins the multi sample variants
# table with a table listing the ACMG genes.

SELECT
refseq_GENE AS Gene,
refseq_CHR AS Chr,
refseq_Tx_START As Transcript_start,
refseq_Tx_END AS Transcript_end,
(refseq_Tx_END - refseq_Tx_START) AS Length_gene,
call.call_set_name AS Sample_id,
COUNT(start) AS Cnt_var
FROM
(
  SELECT
  call.call_set_name,
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
  GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
  GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
  call.FILTER
  FROM
  [_THE_EXPANDED_TABLE_]
  OMIT 
  call IF EVERY(call.FILTER != "PASS")
  OR SOME(call.QC IS NOT NULL)
  HAVING
  cohort_qc IS NULL
) as geno
CROSS JOIN
[_ACMG_GENES_] AS Acmg 
WHERE
geno.reference_name = Acmg.refseq_CHR
AND geno.end <= Acmg.refseq_Tx_END
AND geno.start >= Acmg.refseq_Tx_START
GROUP BY
Gene,
Chr,
Transcript_start,
Transcript_end,
Length_gene,
Sample_id
ORDER BY
Gene,
Sample_id
