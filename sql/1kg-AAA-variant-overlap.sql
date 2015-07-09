SELECT 
AGREEMENT,
AAA,
K1G,
COUNT(*) AS COUNT

FROM (
SELECT 

CASE
  WHEN aaaChr22_allele_frequency >= 0.05 AND k1gChr22_k1000g2012apr_EUR >= 0.05
    THEN "HIGH"
  WHEN aaaChr22_allele_frequency >= 0.005 AND aaaChr22_allele_frequency < 0.05 
    AND k1gChr22_k1000g2012apr_EUR >= 0.005 AND k1gChr22_k1000g2012apr_EUR < 0.05 
    THEN "MODERATE"
  WHEN aaaChr22_allele_frequency >= 0.001 AND aaaChr22_allele_frequency < 0.005 
    AND k1gChr22_k1000g2012apr_EUR >= 0.001 AND k1gChr22_k1000g2012apr_EUR < 0.005 
    THEN "LOW"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR IS NULL
    THEN "NONE"
  ELSE "DISAGREEMENT"
END AS AGREEMENT,

CASE
  WHEN aaaChr22_allele_frequency >= 0.05
    THEN "HIGH"
  WHEN aaaChr22_allele_frequency >= 0.005 AND aaaChr22_allele_frequency < 0.05 
    THEN "MODERATE"
  WHEN aaaChr22_allele_frequency >= 0.001 AND aaaChr22_allele_frequency < 0.005 
    THEN "LOW"  
END AS AAA,

CASE
  WHEN k1gChr22_k1000g2012apr_EUR >= 0.05
    THEN "HIGH"
  WHEN k1gChr22_k1000g2012apr_EUR >= 0.005 AND k1gChr22_k1000g2012apr_EUR < 0.05 
    THEN "MODERATE"
  WHEN k1gChr22_k1000g2012apr_EUR >= 0.001 AND k1gChr22_k1000g2012apr_EUR < 0.005 
    THEN "LOW"  
END AS K1G,

CASE 
  WHEN aaaChr22_allele_frequency IS NOT NULL AND k1gChr22_k1000g2012apr_EUR IS NOT NULL
    THEN "BOTH"
  WHEN aaaChr22_allele_frequency IS NOT NULL AND k1gChr22_k1000g2012apr_EUR IS NULL 
    THEN "AAA"
  WHEN aaaChr22_allele_frequency IS NULL AND k1gChr22_k1000g2012apr_EUR IS NOT NULL
    THEN "1KG"
END AS population,
*
FROM 
  [analysis.aaa_1000g_chr22_AF_full])
WHERE AGREEMENT != "NONE"
#WHERE
#  aaaChr22_allele_frequency IS NOT NULL
#  AND k1gChr22_k1000g2012apr_EUR IS NOT NULL
GROUP BY 
  AGREEMENT,
  AAA,
  K1G
LIMIT 1000