# Determine the value of the Nth quantile for Hardy-Weinberg equilibrium across all variants. 
# 1999th quantile suggested.

SELECT 
  quantile AS cutoff,
  row_num,
FROM(
  SELECT
    quantile,
    ROW_NUMBER() OVER (ORDER BY quantile ASC) row_num,
  FROM (
    SELECT
      QUANTILES(chisq, 2000) AS quantile
    FROM (_MAIN_QUERY_)))
WHERE 
  row_num = _QUANTILE_
