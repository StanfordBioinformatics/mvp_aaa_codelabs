
SELECT
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
  #alternate_allele_count,
  call.DP,
  call.call_set_name,
  FROM(
  
  SELECT
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      #SUM(call.genotype = 1) WITHIN RECORD AS alternate_allele_count,
      call.DP,
      call.call_set_name,
    FROM(
      SELECT
        reference_bases,
        alternate_bases,
        call.genotype,
        call.DP,
        call.call_set_name,
      FROM(FLATTEN((
      [va_aaa_pilot_data.5_genome_test_vcfs_no_calls]), call.DP))
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
      WHERE
        reference_name = 'chr17'
        AND start BETWEEN 41196311
          AND 41277499 )   
    HAVING
    # Skip 1/2 genotypes _and non-SNP variants
    num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
    GROUP BY 
      call.DP,
      call.call_set_name,
    ORDER BY call.DP DESC
