SELECT 
  bc.sample_name AS sample_name,
  (1 - nc.no_call_count/bc.base_count) AS missingness,
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
        [qc_tables.5_genomes_ref_calls_brca1] ) AS ref
    LEFT OUTER JOIN (
      SELECT
        *
      FROM
        [qc_tables.5_genomes_variants_brca1] ) AS var
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
    [qc_tables.5_genomes_ref_calls_brca1]
  GROUP BY
  sample_name ) as bc
ON
  bc.sample_name = nc.sample_name
