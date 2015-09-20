# Compute the ti/tv ratio for known snps grouped by depth of coverage.

SELECT
  call.call_set_name,
  (transitions/transversions) AS titv_ratio,
  average_depth,
  known_snp
FROM (
SELECT
  call.call_set_name,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    ROUND(AVG(call.DP)) AS average_depth,
  IF(dbSNP137 IS NOT NULL, TRUE, FALSE) AS known_snp
FROM (
  SELECT
    reference_name,
    start,
    call.call_set_name,
    mutation,
    call.DP
    
    # Select all calls 
  FROM (
  
    SELECT
      reference_name,
      reference_bases,
      alternate_bases,
      start,
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.DP
    FROM (
      SELECT
        reference_name,
        start,
        call.call_set_name,
        reference_bases,
        alternate_bases,
        call.genotype,
        call.DP,
      FROM(FLATTEN((
        [va_aaa_pilot_data.5_genome_test_vcfs_2]), call.DP))
      # Optionally add clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_  
    )
    WHERE
      call.DP is not null
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')  
  ) 
  GROUP EACH BY
    call.call_set_name,
    reference_name,
    start,
    mutation,
    call.DP
    ) AS all_calls
  
  # Join all the calls against the rsids we found 
  LEFT OUTER JOIN (

      SELECT 
        calls.reference_name AS reference_name,
        start,
        dbSNP137
      FROM (
        SELECT
          reference_name,
          reference_bases,
          alternate_bases,
          start,
          call.call_set_name,
          CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
          COUNT(alternate_bases) WITHIN RECORD AS num_alts,
          call.DP
        FROM (
          SELECT
            reference_name,
            start,
            call.call_set_name,
            reference_bases,
            alternate_bases,
            call.genotype,
            call.DP,
          FROM ( 
            FLATTEN ((
              [va_aaa_pilot_data.5_genome_test_vcfs_2]), call.DP
            )
          )
          # Optionally add clause here to limit the query to a particular
          # region of the genome.
          #_WHERE_  
        )
        WHERE
          call.DP is not null
        HAVING
          # Skip 1/2 genotypes _and non-SNP variants
          num_alts = 1
          AND reference_bases IN ('A','C','G','T')
          AND alternate_bases IN ('A','C','G','T')) AS calls
        JOIN EACH 
        (SELECT 
          CONCAT('chr', STRING(CHROM)) as chr,
          POS,
          dbSNP137,
        FROM 
          [1KGenomeAnnotation.annotations] 
        WHERE
          dbSNP137 IS NOT NULL) as annotation
      ON
        calls.reference_name = annotation.chr
        AND calls.start = annotation.POS
      GROUP BY 
        reference_name,
        start,
        dbSNP137 
    ) AS rsids 
        
    ON 
      all_calls.reference_name = rsids.reference_name
      AND all_calls.start = rsids.start 

GROUP EACH BY
  call.call_set_name,
  known_snp,
  call.DP
)
WHERE
  transversions > 0
GROUP BY
  call.call_set_name,
  average_depth,
  known_snp,
  titv_ratio,
  call.DP
ORDER BY average_depth DESC
 LIMIT 1000
