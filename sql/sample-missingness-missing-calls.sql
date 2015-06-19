SELECT 
  bc.sample_name AS sample_name,
  nc.no_call_count as no_call_count,
  nc.no_call_count/bc.base_count AS missingness,
FROM (
  SELECT 
    sample_name,
    COUNT(no_call) as no_call_count,
  FROM(
    SELECT
      ref.sample_name AS sample_name,
      ref.chromosome AS chromosome,
      ref.start AS start,
      ref.end AS end,
      ref.reference_bases AS reference_bases,
      True AS no_call
    FROM
      (SELECT 
        *
       FROM
        [_REF_TABLE_] ) AS ref
    LEFT OUTER JOIN (
      SELECT
        *
      FROM
        [_VARIANT_TABLE_] ) AS var
    ON
      ref.sample_name = var.sample_name
      AND ref.chromosome = var.chromosome
      AND ref.start = var.start
    WHERE
      ref.is_ref is False
      AND var.is_variant_call is NULL )
  GROUP BY
    sample_name ) AS nc
JOIN (
  SELECT
    sample_name,
    COUNT(start) AS base_count
  FROM 
    [_REF_TABLE_]
  GROUP BY
  sample_name ) as bc
ON
  bc.sample_name = nc.sample_name;
