SELECT 
  sample_id,
  gender,
  perct_het_alt_in_snvs,
FROM (
SELECT 
  data.sample_id AS sample_id,
  data.perct_het_alt_in_snvs AS perct_het_alt_in_snvs,
  data.all_callable_sites,
  data.hom_AA_count,
  data.het_RA_count,
  data.hom_RR_count,
  data.all_snvs,
  info.SEX AS gender,
FROM (
SELECT
  call.call_set_name AS sample_id,
  ROUND((het_RA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_het_alt_in_snvs,
  ROUND((hom_AA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_hom_alt_in_snvs,
  (hom_AA_count + het_RA_count + hom_RR_count) AS all_callable_sites,
  hom_AA_count,
  het_RA_count,
  hom_RR_count,
  (hom_AA_count + het_RA_count) AS all_snvs,
FROM
  (
  SELECT
    call.call_set_name,
    SUM(0 = first_allele
      AND 0 = second_allele) AS hom_RR_count,
    SUM(first_allele = second_allele AND first_allele > 0) AS hom_AA_count,
    SUM((first_allele != second_allele OR second_allele IS NULL)
      AND (first_allele > 0 OR second_allele > 0)) AS het_RA_count
  FROM (
    SELECT
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
    FROM
      [_THE_TABLE_]
    WHERE
      reference_name = 'chrX'
      AND start NOT BETWEEN 59999 AND 2699519
      AND start NOT BETWEEN 154931042 AND 155260559
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    call.call_set_name)) AS data
JOIN (
  SELECT
    IlluminaID,
    SEX
  FROM
    [_PATIENT_INFO_]) AS info
ON
  data.sample_id = info.IlluminaID)
WHERE
  (gender = 'M' AND perct_het_alt_in_snvs > _MALE_CUTOFF_) OR
  (gender = 'F' AND perct_het_alt_in_snvs < _FEMALE_CUTOFF_) 
ORDER BY
  data.sample_id
