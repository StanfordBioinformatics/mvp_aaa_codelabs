# Identify heterozygous variants on the sex chromosomes for reportedly male samples. 

SELECT
  variant_id,
  sample_id,
  "heterozygous_haplotype" AS failure_reason,
FROM (
SELECT
  variant_id,
  reference_name,
  start,
  end,
  reference_bases,
  call.call_set_name AS sample_id,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
FROM(FLATTEN((
  [_THE_EXPANDED_TABLE_]), call.call_set_name))
WHERE
  reference_name IN ('chrX', 'chrY')
OMIT
  call IF (2 > COUNT(call.genotype))
  OR EVERY(call.genotype <= 0)
  OR EVERY(call.genotype = 1)
  # Pseudoautosomal Region 1
  OR (reference_name = 'chrX'
    AND start > 60001
    AND end < 2699520)
  OR (reference_name = 'chrY'
    AND start > 10001
    AND end < 2649520)
  # Pseudoautosomal Region 2
  OR (reference_name = 'chrX'
    AND start > 155260560
    AND end < 155270560)
  OR (reference_name = 'chrY' 
    AND start > 59363566
    AND end < 59373566)) AS seq
JOIN (
  SELECT
    IlluminaID,
    SEX
  FROM
    [_PATIENT_INFO_] ) AS info
ON
  seq.sample_id = info.IlluminaID
WHERE
  SEX = 'M'
# Optionally add a clause here to sort and limit the results.
#_ORDER_BY_
