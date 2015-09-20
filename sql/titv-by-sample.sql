# Calculate the ti/tv ratio for each sample.

SELECT
call.call_set_name AS sample.id,
transitions,
transversions,
transitions/transversions AS titv,
num_variants,
FROM 
(
  SELECT
  call.call_set_name,
  SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
  SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                   'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
  COUNT(mutation) AS num_variants
  FROM 
  (
    SELECT
    call.call_set_name,
    CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    GROUP_CONCAT(QC) WITHIN RECORD AS qc,
    GROUP_CONCAT(call.QC) WITHIN CALL AS call_qc,
    FROM
    FLATTEN([_THE_EXPANDED_TABLE_], alternate_bases)
    OMIT call IF EVERY (call.genotype <= 0)
    OR SOME(call.qc IS NOT NULL)
    HAVING
    num_alts = 1
    AND QC IS NULL
  )
  GROUP BY
  call.call_set_name
)
