SELECT 
CHROM,
POS,
ID,
REF,
ALT,
QUAL,
FILTER,
INFO,
sample_calls
FROM (
  SELECT
  CHROM,
  POS,
  ID,
  REF,
  ALT,
  ROUND(AVG(QUAL),2) AS QUAL,
  FILTER,
  CONCAT("DP=", STRING(ROUND(AVG(DP),2))) AS INFO,
  COUNT(sample) AS sample_count,
  GROUP_CONCAT(sample) AS sample_calls
  FROM (
    SELECT
    reference_name AS CHROM,
    start AS POS,
    names AS ID,
    reference_bases AS REF,
    alternate_bases AS ALT,
    AVG(quality) AS QUAL,
    filter AS FILTER,
    DP,
    CONCAT(call_set_name, ":", genotype) AS sample
    FROM (
      SELECT
      reference_name,
      start,
      GROUP_CONCAT(names) WITHIN RECORD AS names,
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      quality,
      GROUP_CONCAT(filter) WITHIN RECORD AS filter, 
      call.DP AS DP,
      call.call_set_name AS call_set_name,
      GROUP_CONCAT(STRING(call.genotype),"/") WITHIN call AS genotype,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      FROM
      [_THE_EXPANDED_TABLE_] 
      WHERE
      reference_name = _CHR_
      OMIT 
      call IF EVERY(call.genotype < 0) 
      HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      AND names is not null
      AND FILTER == 'PASS'
    )
    GROUP EACH BY 
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    FILTER,
    sample,
    DP)
  GROUP EACH BY 
  CHROM,
  POS,
  ID,
  REF,
  ALT,
  FILTER,)
# Only include sites where all samples have a call
WHERE
sample_count = 479
#ORDER BY
# POS,
# ALT