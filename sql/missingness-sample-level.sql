# Determine the level of missingness for each sample when compared to the hg19 reference genome
SELECT
  sample_id,
  ROUND(1 - ((build_length - all_calls_count)/build_length), 3) AS missingness,
FROM (
  SELECT
    call.call_set_name AS sample_id,
    ref_count + alt_count AS all_calls_count,
    ref_count,
    alt_count,
    3088269832 AS build_length,
  FROM (
    SELECT
      call.call_set_name,
      SUM(IF(genotype = '0,0',
        (end - start),
        0)) AS ref_count,
      SUM(IF(genotype NOT IN ('0,0', '-1,-1'),
        1,
        0)) AS alt_count,
    FROM (
    SELECT
      call.call_set_name,
      start,
      end,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alts,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
    FROM
      [gbsc-gcp-project-mvp:va_aaa_pilot_data.aaa_20160709_genome_calls_no_qc]
    WHERE
      reference_bases IN ('A', 'C', 'T', 'G')
    HAVING
      alts IN ('A', 'C', 'T', 'G', '<NON_REF>'))
    GROUP BY
      call.call_set_name)) 
ORDER BY
  missingness DESC
