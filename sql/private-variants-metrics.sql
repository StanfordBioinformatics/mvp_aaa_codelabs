# Compute the average and standard deviation of private variants within the cohort.
# This query requires that _main_query_ be substituted with the contents of private-snv-counts.sql.
SELECT 
ROUND(AVG(private_variant_count), 3) AS average,
ROUND(STDDEV(private_variant_count), 3) AS stddev
FROM (_MAIN_QUERY_)
