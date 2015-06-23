# Compute private variants counts for each sample.
SELECT 
ROUND(AVG(private_variant_count), 3) AS average,
ROUND(STDDEV(private_variant_count), 3) AS stddev
FROM (_MAIN_QUERY_)))
