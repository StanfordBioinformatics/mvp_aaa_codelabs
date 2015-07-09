SELECT
  frequency,
  population,
  COUNT(CONCAT(population, frequency)) AS count
FROM (
SELECT 
*,
CASE
  WHEN aaaChr22_allele_frequency >= 0.05 AND k1gChr22_k1000g2012apr_EUR >= 0.05
    THEN "HIGH"
  WHEN aaaChr22_allele_frequency >= 0.05 AND k1gChr22_k1000g2012apr_EUR IS NULL
    THEN "HIGH"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR >= 0.05
    THEN "HIGH"
  WHEN aaaChr22_allele_frequency BETWEEN 0.005 AND 0.05 AND k1gChr22_k1000g2012apr_EUR BETWEEN 0.005 AND 0.05 
    THEN "MODERATE"
  WHEN aaaChr22_allele_frequency BETWEEN 0.005 AND 0.05 AND k1gChr22_k1000g2012apr_EUR IS NULL
    THEN "MODERATE"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR BETWEEN 0.005 AND 0.05 
    THEN "MODERATE" 
  WHEN aaaChr22_allele_frequency BETWEEN 0.001 AND 0.005 AND k1gChr22_k1000g2012apr_EUR BETWEEN 0.001 AND 0.005 
    THEN "LOW"
  WHEN aaaChr22_allele_frequency BETWEEN 0.001 AND 0.005 AND k1gChr22_k1000g2012apr_EUR IS NULL
    THEN "LOW"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR BETWEEN 0.001 AND 0.005 
    THEN "LOW"    
END AS frequency,

CASE 
  WHEN aaaChr22_allele_frequency IS NOT NULL AND k1gChr22_k1000g2012apr_EUR IS NOT NULL
    THEN "BOTH"
  WHEN aaaChr22_allele_frequency IS NOT NULL AND k1gChr22_k1000g2012apr_EUR IS NULL 
    THEN "AAA"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR IS NOT NULL
    THEN "1KG"
END AS population,

FROM 
  [analysis.aaa_1000g_chr22_AF_comparison_20150630_pos_corrected] )
WHERE
  population IS NOT NULL 
  AND frequency IS NOT NULL
GROUP BY 
  frequency,
  population
