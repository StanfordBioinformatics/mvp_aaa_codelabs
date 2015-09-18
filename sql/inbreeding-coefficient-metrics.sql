# Calculate the average and standard deviation for the inbreeding coefficient.
# This query requires that _MAIN_QUERY_ be substituted with the contents of homozygous-variants.sql

SELECT
ROUND(AVG(F), 3) AS average,
ROUND(STDDEV(F), 3) AS stddev
FROM (_MAIN_QUERY_)
