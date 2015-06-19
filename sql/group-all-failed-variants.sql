SELECT 
variant_id,
failure_reason,
sample_id
FROM (
  SELECT
  variant_id,
  failure_reason,
  sample_id,
  FROM
  [qc_tables.heterozygous_haplotype],
  [qc_tables.titv_depth]),
(SELECT
 variant_id,
 failure_reason,
 "" AS sample_id
 FROM
 [qc_tables.blacklisted], 
 [qc_tables.hardy_weinberg], 
 [qc_tables.titv_genomic_window], 
 [qc_tables.variant_missingness])