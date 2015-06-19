SELECT
variant_id,
reference_name,
start,
depth
FROM(
  SELECT
  variant_id,
  reference_name,
  start,
  ROUND(AVG(call.DP)) AS depth
  FROM
  [_THE_EXPANDED_TABLE_]
  GROUP EACH BY 
  variant_id,
  reference_name,
  start)
WHERE
depth < _MIN_ OR
depth > _MAX_
