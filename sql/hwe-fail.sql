# Get all variants that have a chi squared value for Hardy-Weinberg equilibrium above a definited limit
SELECT
  variant_id,
  chisq,
  "hardy_weinberg" AS failure_reason,
FROM (_MAIN_QUERY_)
WHERE
  chisq > _CUTOFF_
#_ORDER_BY_

