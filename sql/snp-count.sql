# overly complicated but does the job

SELECT
  COUNT(mutation) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([va_aaa_pilot_data.5_genome_test_gvcfs_2], call.call_set_name)
      OMIT RECORD if every(call.genotype <= 0)
    HAVING
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
    )
  WHERE
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype)