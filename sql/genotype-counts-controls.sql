# Count the number of homozygous negative, heterozygous, homozygous positive, and multiallelic positions
# in genomes NOT associated with AAA (aaa size < 3.0).  Output format is based on requirements from
# recipient (Dana).

SELECT
chromosome,
start,
# Genotype counts
SUM(CASE WHEN genotype = "[0,1]" AND aaa_size < 3 THEN 1 ELSE 0 END) AS nr_sample_het01_SNV,
SUM(CASE WHEN genotype = "[1,2]" AND aaa_size < 3 THEN 1 ELSE 0 END) AS nr_sample_het12_SNV,
SUM(CASE WHEN genotype = "[1,1]" AND aaa_size < 3 THEN 1 ELSE 0 END) AS nr_sample_hom11_SNV,
SUM(CASE WHEN genotype = "[0,0]" AND aaa_size < 3 THEN 1 ELSE 0 END) AS nr_sample_hom00_REF,
FROM (
  
  SELECT
  seq.sample_id AS sample_id,
  seq.reference_name AS chromosome,
  seq.start AS start,
  seq.end AS end,
  seq.genotype AS genotype,
  phen.AAA_SIZE AS aaa_size,
  FROM (
    SELECT
    sample_id,
    reference_name,
    start,
    end,
    genotype,
    FROM js(
      (SELECT
       call.call_set_name,
       reference_name,
       start,
       end,
       call.genotype,
       reference_bases,
       GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alts,
       COUNT(alternate_bases) WITHIN RECORD AS num_alts,
       FROM
       [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
       WHERE
       reference_name = 'chr17'
       AND start BETWEEN 41196311
       AND 41277499
       OMIT 
       call IF EVERY (call.genotype < 0)
       HAVING 
       reference_bases IN ('A','C','G','T')
      ),
      // Start javascript function
      // Input Columns
      call.call_set_name, reference_name, start, end, call.genotype,
      // Output Schema
      "[{name: 'sample_id', type: 'string'},
        {name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'end', type: 'integer'},
        {name: 'genotype', type: 'string'}]",
      // Function
      "function(r, emit) {
          for (c of r.call) {
            var genotype = JSON.stringify(c.genotype.sort());
              emit({
                sample_id: c.call_set_name,
                reference_name: r.reference_name,
                start: r.start,
                end: r.end,
                genotype: genotype,
              })
          }
        }")) AS seq
  JOIN EACH (
    SELECT
    IlluminaID,
    AAA_SIZE,
    FROM
    va_aaa_pilot_data.patient_info ) AS phen
  ON
  seq.sample_id = phen.IlluminaID)
GROUP BY
chromosome,
start,

ORDER BY 
chromosome,
start
