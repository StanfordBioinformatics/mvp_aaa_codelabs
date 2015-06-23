# Compute private variants counts for each sample.
SELECT 
ROUND(AVG(private_variant_count), 3) AS average,
ROUND(STDDEV(private_variant_count), 3) AS stddev
FROM (
  SELECT
  call.call_set_name,
  COUNT(call.call_set_name) AS private_variant_count,
  FROM (
    SELECT
    reference_name,
    start,
    GROUP_CONCAT(CASE WHEN cnt = 1 THEN 'S'
                 WHEN cnt = 2 THEN 'D'
                 ELSE STRING(cnt) END) AS SINGLETON_DOUBLETON,
    reference_bases,
    alternate_bases,
    GROUP_CONCAT(call.call_set_name) AS call.call_set_name,
    GROUP_CONCAT(genotype) AS genotype,
    SUM(num_samples_with_variant) AS num_samples_with_variant
    FROM (_MAIN_QUERY_)))
