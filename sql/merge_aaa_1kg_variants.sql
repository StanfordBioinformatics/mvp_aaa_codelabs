# Merge variants from the AAA dataset and the 1000 Genomes dataset.
# NOTE: some counts are hard-coded in this query.  Use with caution.

SELECT *
  FROM
(
  SELECT 
  reference_name, 
  start,
  end,
  reference_bases,
  alternates,
  num_samples,
  ROUND(num_samples/446, 4) AS allele_frequency,
  INTEGER(FLOOR(start / 100000)) AS window,
  FROM (
    SELECT
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
    FROM
    [gbsc-gcp-project-mvp:va_aaa_pilot_data.multi_sample_variants_full_qc]
    OMIT 
    call IF EVERY(call.FILTER != "PASS")
    OR SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL 
  )
  GROUP EACH BY
  reference_name,
  start,
  end,
  reference_bases,
  alternates,
  call_qc,
  cohort_qc,
  num_samples,
  allele_frequency,
  window
) AS aaaChr22
FULL OUTER JOIN EACH
(
  SELECT
  CONCAT(STRING("chr"), STRING(CHROM)) AS Chr,
  POS-1 AS Start,
  REF,
  ALT,
  k1000g2012apr_EUR,
  INTEGER(FLOOR((POS-1) / 100000)) AS window,
  FROM 
  [stanford.edu:gbsc-stanford-google:1KGenomeAnnotation.annotations] 
) AS k1gChr22
ON
aaaChr22.reference_name = k1gChr22.Chr
AND aaaChr22.window = k1gChr22.window
AND aaaChr22.alternates = k1gChr22.ALT

WHERE
aaaChr22.start = k1gChr22.Start

