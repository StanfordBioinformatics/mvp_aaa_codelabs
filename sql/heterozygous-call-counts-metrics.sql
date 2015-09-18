# Compute the average and standard deviation of heterzygous call counts for 
# all genomes.  This query requires that _MAIN_QUERY_ be substituted with
# the contents of heterozygous-call-counts.sql

SELECT
ROUND(AVG(O_HET), 3) AS average,
ROUND(STDDEV(O_HET), 3) AS stddev
FROM (_MAIN_QUERY_)
