 SELECT
  vars.call_set_name,
  ref.chromosome,
  ref.start,
  ref.end,
  ref.reference_bases,
  vars.var_start,
  vars.var_end,
  
  #vars.variant_called_count,
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
        #SUM(called_count) AS variant_called_count,
      FROM ( FLATTEN((
        # _LIMIT the query to SNPs
        SELECT
          call.call_set_name,
          reference_name as chromosome,
          start,
          end,
          reference_bases,
          LENGTH(reference_bases) AS ref_len,
          MIN(LENGTH(alternate_bases)) WITHIN RECORD AS alt_len,
          IF(alternate_bases IS NULL,
            FALSE,
            TRUE) AS is_variant_call,
          #SUM(call.genotype >= 0) WITHIN RECORD AS called_count,
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
        reference_bases ) AS vars
ON 
  ref.chromosome = vars.chromosome
  AND ref.start = vars.var_start #) AS
  
ORDER BY 
  start,
  end,
LIMIT 1000        
