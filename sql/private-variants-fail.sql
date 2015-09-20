# Identify samples whose private variant count is outside of a defined range.
# This query requires that _MAIN_QUERY_ be substituted with the contents of private-snv-counts.sql.

SELECT 
call.call_set_name AS sample_id,
private_variant_count
FROM (_MAIN_QUERY_)
WHERE
private_variant_count  > _MAX_VALUE_ OR
private_variant_count < _MIN_VALUE_
