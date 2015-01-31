# Compute the Ti/Tv ratio of the 1,000 Genomes dataset.
SELECT
  transitions,
  transversions,
  transitions/transversions AS titv,
  COUNT
FROM (
  SELECT
    SUM(IF(mutation IN ('A->G',
          'G->A',
          'C->T',
          'T->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transitions,
    SUM(IF(mutation IN ('A->C',
          'C->A',
          'G->T',
          'T->G',
          'A->T',
          'T->A',
          'C->G',
          'G->C'),
        INTEGER(num_snps),
        INTEGER(0))) AS transversions,
        COUNT(mutation) AS COUNT
  FROM (
    SELECT
      CONCAT(reference_bases,
        CONCAT(STRING('->'),
          alternate_bases)) AS mutation,
      COUNT(alternate_bases) AS num_snps,
    FROM
      FLATTEN([_THE_TABLE_], call)
    WHERE
      LENGTH(alternate_bases) == 1
      AND LENGTH(reference_bases) == 1
      _LIMITS_
    GROUP BY
      mutation,
    ORDER BY
      mutation))