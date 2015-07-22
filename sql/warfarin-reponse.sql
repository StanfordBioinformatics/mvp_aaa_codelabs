SELECT
call.call_set_name,
names,
reference_name,
start,
end,
genotype,
FROM (
  SELECT
  call.call_set_name,
  reference_name,
  start,
  end,
  names,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
  FROM
  (FLATTEN((va_aaa_pilot_data.multi_sample_variants_full_qc),call.call_set_name))
  HAVING
  names in ('rs1057910', 'rs1799853', 'rs9923231')
)
GROUP BY
call.call_set_name,
names,
reference_name,
start,
end,
names,
genotype