SELECT
 refs.sample_name AS sample_name,
 refs.chromosome AS chromosome,
 refs.start AS start,
 refs.end AS end,
 refs.reference_bases AS reference_bases,
 IF(calls.ref_block_start IS NULL, False, True) AS is_ref,
FROM (
  SELECT
    sample_names.call_set_name AS sample_name,
    ref.chromosome AS chromosome,
    ref.start AS start,
    ref.END AS end,
    ref.reference_bases AS reference_bases,
  FROM (
    SELECT
      call.call_set_name AS call_set_name,
    FROM
      [va_aaa_pilot_data.5_genome_test_gvcfs]
    GROUP BY
      call_set_name
  ) as sample_names

  CROSS JOIN (
  
    SELECT  
      CONCAT('chr', Chr) AS chromosome,
      Start,
      END,
      Ref AS reference_bases,
    FROM 
      [google.com:biggene:test.hg19] 
    WHERE 
      Chr = '17'
      AND Start BETWEEN 41196311
        AND 41277499 ) AS ref ) AS refs 

LEFT OUTER JOIN (
 
 SELECT
  refs.sample_name AS sample_name,
  refs.chromosome AS chromosome,
  refs.start AS start,
  refs.end AS end,
  refs.reference_bases AS reference_bases,
  ref_calls.ref_block_start AS ref_block_start,
  ref_calls.ref_block_end AS ref_block_end,
FROM (
  SELECT
    sample_names.call_set_name AS sample_name,
    ref.chromosome AS chromosome,
    ref.start AS start,
    ref.END AS end,
    ref.reference_bases AS reference_bases,
  FROM (
    SELECT
      call.call_set_name AS call_set_name,
    FROM
      [va_aaa_pilot_data.5_genome_test_gvcfs]
    GROUP BY
      call_set_name
  ) as sample_names

  CROSS JOIN (
  
    SELECT  
      CONCAT('chr', Chr) AS chromosome,
      Start,
      END,
      Ref AS reference_bases,
    FROM 
      [google.com:biggene:test.hg19] 
    WHERE 
      Chr = '17'
      AND Start BETWEEN 41196311
        AND 41277499 ) AS ref ) AS refs 

    LEFT OUTER JOIN (
    
      SELECT
        call.call_set_name AS call_set_name,
        reference_name AS chromosome,
        start AS ref_block_start,
        end AS ref_block_end,
      IF(alternate_bases IS NULL,
        FALSE,
        TRUE) AS is_variant_call,
      FROM( FLATTEN((
        SELECT
          call.call_set_name,
          reference_name,
          start,
          end,
          alternate_bases,
        FROM  
          [va_aaa_pilot_data.5_genome_test_gvcfs]), call.call_set_name))
        WHERE 
          reference_name = 'chr17'
          AND start BETWEEN 41196311
          AND 41277499 ) AS ref_calls
    ON 
      refs.chromosome = ref_calls.chromosome
      AND refs.sample_name = ref_calls.call_set_name
    WHERE
      refs.start >= ref_calls.ref_block_start 
      AND refs.end <= ref_calls.ref_block_end ) AS calls
ON
  refs.chromosome = calls.chromosome
  AND refs.sample_name = calls.sample_name
  AND refs.start = calls.start
ORDER BY 
  refs.start,
  refs.end,
  refs.sample_name,
     
