# Count the number of indels by length
SELECT
  length_alt,
  length_ref,
  COUNT(length_alt) AS count,
FROM(
  SELECT
    mutation,
    COUNT(mutation) AS num_variants,
    start,
    end,
    genotype,
    length_alt,
    length_ref,
  FROM (
    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      LENGTH(alternate_bases) AS length_alt,
      LENGTH(reference_bases) AS length_ref,
      reference_name,
      start,
      end,
      GROUP_CONCAT(STRING(call.genotype)) WITHIN RECORD AS genotype,
    FROM
      FLATTEN([_THE_TABLE_], call.call_set_name)
    OMIT RECORD IF EVERY (call.genotype <= 0)
    HAVING
      num_alts = 1
      AND LENGTH(reference_bases) > 1
      OR LENGTH(alternate_bases) > 1
  )
  WHERE 
    genotype != '-1,-1'
    #_AND_
  GROUP EACH BY
    mutation,
    start,
    end,
    genotype,
    length_ref,
    length_alt)
GROUP EACH BY 
  length_alt,
  length_ref,
ORDER BY 
  length_ref,
  length_alt DESC
