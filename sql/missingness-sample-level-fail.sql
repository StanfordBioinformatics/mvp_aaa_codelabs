# Identify genomes where the missingness value is above a defined cutoff (0.9 recommended).

SELECT
  sample_id,
  missingness
FROM (
SELECT
  sample_id,
  ROUND(1 - ((build_length - all_calls_count)/build_length), 3) AS missingness,

FROM (
  SELECT
    call.call_set_name AS sample_id,
    ref_count + alt_count AS all_calls_count,
    ref_count,
    alt_count,
    3049315783 AS build_length,
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
      [_THE_TABLE_]
    WHERE
      reference_bases IN ('A', 'C', 'T', 'G')
    HAVING
      alts IN ('A', 'C', 'T', 'G', '<NON_REF>'))
    GROUP BY
      call.call_set_name)))
WHERE
  missingness < _CUTOFF_
#ORDER BY
#  g.sample_id
