# Select genomes that are outside a defined range for the number of heterozygous calls.  
# This query requires that _main_query_ be substituted with the contents of heterozygous-call-counts.sql

SELECT
call.call_set_name AS sample_id,
O_HET
FROM (_MAIN_QUERY_)
WHERE
O_HET > _MAX_VALUE_ OR 
O_HET < _MIN_VALUE_
