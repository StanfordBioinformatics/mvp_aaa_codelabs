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
