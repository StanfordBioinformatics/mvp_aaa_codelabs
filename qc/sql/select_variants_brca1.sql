SELECT
  vars.call_set_name AS sample_name,
  ref.chromosome AS chromosome,
  ref.start AS start,
  vars.var_end AS end,
  ref.reference_bases AS reference_bases,
  vars.alternate_bases AS alternate_bases,
  vars.is_variant_call AS is_variant_call,
FROM (
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
      AND 41277499 ) AS ref 
JOIN (
SELECT
        call.call_set_name AS call_set_name,
        chromosome,
        start AS var_start,
        end AS var_end,
        INTEGER(FLOOR(start / 5000)) AS bin,
        reference_bases,
        alternate_bases,
        is_variant_call,
      FROM ( FLATTEN((
        SELECT
          call.call_set_name,
          reference_name as chromosome,
          start,
          end,
          reference_bases,
          LENGTH(reference_bases) AS ref_len,
          MIN(LENGTH(alternate_bases)) WITHIN RECORD AS alt_len,
          alternate_bases,
          IF(alternate_bases IS NULL,
            FALSE,
            TRUE) AS is_variant_call,
        FROM
          [va_aaa_pilot_data.5_genome_test_gvcfs]
        WHERE
          reference_name = 'chr17'
          AND start BETWEEN 41196311
            AND 41277499                
        HAVING
          ref_len = 1
          AND alt_len = 1
          AND is_variant_call), call.call_set_name))
      GROUP EACH BY
        call_set_name,
        chromosome,
        var_start,
        var_end,
        bin,
        reference_bases,
        alternate_bases,
        is_variant_call ) AS vars
ON 
  ref.chromosome = vars.chromosome
  AND ref.start = vars.var_start
ORDER BY 
  start,
  end  
