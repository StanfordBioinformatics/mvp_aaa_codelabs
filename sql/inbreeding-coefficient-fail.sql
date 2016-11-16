# Identify genomes that have an inbreeding coefficient outside a defined range.
# This query requires that _main_query_ be replaced with the contents of homozygous-variants.sql

SELECT
call.call_set_name AS sample_id,
O_HOM,
E_HOM,
N_SITES,
F
FROM (_MAIN_QUERY_)
WHERE
F > _MAX_VALUE_ OR 
F < _MIN_VALUE_
