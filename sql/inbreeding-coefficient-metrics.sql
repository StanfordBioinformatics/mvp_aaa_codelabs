# Compute the expected and observed homozygosity rate for each individual.
SELECT
ROUND(AVG(F), 3) AS average,
ROUND(STDDEV(F), 3) AS stddev
FROM (_MAIN_QUERY_)
