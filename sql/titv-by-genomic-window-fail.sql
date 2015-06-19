SELECT
var.variant_id AS variant_id,
titv,
"titv_by_genomic_window" AS failure_reason,
FROM (
  SELECT
  variant_id,
  reference_name,
  start,
  end,
  INTEGER(FLOOR(start / 100000)) AS window,
  FROM
  [_THE_EXPANDED_TABLE_]
  OMIT call IF EVERY(call.genotype <= 0)
  # Optionally add clause here to limit the query to a particular
  # region of the genome.
  #_WHERE_ 
) AS var
JOIN (
  SELECT
  reference_name,
  window,
  window_start,
  transitions,
  transversions,
  titv,
  num_variants_in_window,
  FROM (
    SELECT
    reference_name,
    window,
    window * 100000 AS window_start,
    transitions,
    transversions,
    transitions/transversions AS titv,
    num_variants_in_window,
    FROM (
      SELECT
      reference_name,
      window,
      SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
      SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                       'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
      COUNT(mutation) AS num_variants_in_window
      FROM (
        SELECT
        reference_name,
        INTEGER(FLOOR(start / 100000)) AS window,
        CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        FROM
        [_THE_EXPANDED_TABLE_]
        # Optionally add clause here to limit the query to a particular
        # region of the genome.
        WHERE reference_name = 'chr22'
        HAVING
        # Skip 1/2 genotypes _and non-SNP variants
        num_alts = 1
        AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T'))
      GROUP BY
      reference_name,
      window))
  WHERE
  titv > _MAX_ OR
  titv < _MIN_) as win
ON 
var.window = win.window
GROUP BY 
variant_id,
var.reference_name,
win.window_start,
titv,
#_ORDER_BY_

