# Identify the alternate allele counts for which the Ti/Tv ratio is outside a defined range.
SELECT
  *
FROM (
  SELECT
    transitions,
    transversions,
    transitions/transversions AS titv,
    alternate_allele_count
  FROM (
    SELECT
      SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
      SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                       'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
      alternate_allele_count
      FROM (
        SELECT
          CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
          COUNT(alternate_bases) WITHIN RECORD AS num_alts,
          SUM(call.genotype = 1) WITHIN RECORD AS alternate_allele_count,
        FROM
          [_THE_EXPANDED_TABLE_]
        # Optionally add clause here to limit the query to a particular
        # region of the genome.
        #_WHERE_
        OMIT 
          call IF EVERY(call.genotype <= 0)
        HAVING
        # Skip 1/2 genotypes _and non-SNP variants
          num_alts = 1
          AND reference_bases IN ('A','C','G','T')
          AND alternate_bases IN ('A','C','G','T'))
      GROUP BY
        alternate_allele_count))
WHERE
  titv > _MAX_ OR
  titv < _MIN_
#_ORDER_BY_
