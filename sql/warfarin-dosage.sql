# Determine the number of samples that require each Warfarin dosage recommmendation based on their genotypes.

SELECT
  dosage,
  COUNT(*) AS count
FROM (
SELECT 
  call.call_set_name,
  snps,
  CASE
    WHEN (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:ref%') THEN '5-7'
    WHEN (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:ref%') THEN '5-7' 
    WHEN (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:ref%') THEN '3-4'
    
    WHEN (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1057910:ref%' AND snps LIKE '%rs1799853:het%') THEN '5-7'
    WHEN (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1057910:ref%' AND snps LIKE '%rs1799853:het%') THEN '3-4'
    WHEN (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1057910:ref%' AND snps LIKE '%rs1799853:het%') THEN '3-4' 
    
    WHEN (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:het%') THEN '3-4'
    WHEN (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:het%') THEN '3-4' 
    WHEN (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:het%') THEN '0.5-2' 

    WHEN (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:ref%') THEN '3-4'
    WHEN (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:ref%') THEN '3-4' 
    WHEN (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:ref%') THEN '0.5-2' 
    
    WHEN ((snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:het%' AND snps LIKE '%rs1057910:het%')
      OR (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:hom%')) THEN '3-4'
    WHEN ((snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:het%' AND snps LIKE '%rs1057910:het%')
      OR (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:hom%')) THEN '0.5-2'
    WHEN ((snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:het%' AND snps LIKE '%rs1057910:het%')
      OR (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:hom%' AND snps LIKE '%rs1057910:hom%')) THEN '0.5-2'
    
    WHEN (snps LIKE '%rs9923231:ref%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:hom%') THEN '0.5-2'
    WHEN (snps LIKE '%rs9923231:het%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:hom%') THEN '0.5-2' 
    WHEN (snps LIKE '%rs9923231:hom%' AND snps LIKE '%rs1799853:ref%' AND snps LIKE '%rs1057910:hom%') THEN '0.5-2'
    ELSE 'undeterminable'
  END AS dosage
FROM (
SELECT
  call.call_set_name,
  GROUP_CONCAT(geno) AS snps
FROM (
SELECT
  call.call_set_name,
  names,
  genotype,
  CASE 
    WHEN (genotype = '0,0') THEN CONCAT(names, ':ref')
    WHEN (genotype = '0,1') THEN CONCAT(names, ':het')
    WHEN (genotype = '1,1') THEN CONCAT(names, ':hom')
    ELSE NULL
  END AS geno
FROM (
SELECT
  call.call_set_name,
  names,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
FROM
  (FLATTEN((_THE_EXPANDED_TABLE_),call.call_set_name))
HAVING
  names in ('rs1057910', 'rs1799853', 'rs9923231')
)
GROUP BY
  call.call_set_name,
  names,
  genotype,
  geno)
GROUP BY
  call.call_set_name))
GROUP BY
  dosage
