SELECT
call.call_set_name,
reference_name,
contig.len,
num_SNVs,
num_REFs,
(num_SNVs + num_REFs) AS num_called_point_pos,
ROUND((num_SNVs + num_REFs) / contig.len, 3) AS prop_w_point_ino,
(contig.len - (num_SNVs + num_REFs)) AS pos_no_point_info,
ROUND((contig.len - (num_SNVs + num_REFs)) / contig.len, 3) prop_no_point_info
FROM
(
  SELECT
  call.call_set_name,
  reference_name,
  assembly.LENGTH AS contig.len,
  SUM(call.FILTER="PASS" AND (LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")))) AS num_SNVs,
  SUM(IF (genotypes=="0/0", (end - start), 0)) AS num_REFs
  FROM (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotypes,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    call.FILTER
    FROM 
    [_THE_TABLE_]
    OMIT
    call IF SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL
  ) AS geno
  JOIN 
  [stanford.edu:gbsc-stanford-google:resources.hg19_Assembly_BinaRuns] AS assembly
  ON
  geno.reference_name = assembly.CHR
  GROUP BY
  call.call_set_name,
  reference_name,
  contig.len
)
ORDER BY
call.call_set_name,
reference_name,