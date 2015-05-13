SELECT
  seq.variant_id AS variant_id,
  seq.reference_name AS reference_name,
  seq.start AS start,
  seq.end AS end,
  bl.Artifact_Type AS Artifact_Type 
FROM (
  SELECT
    variant_id,
    reference_name,
    start,
    end,
  FROM
    [_THE_EXPANDED_TABLE_]) as seq
JOIN (
  SELECT 
    reference_name,
    start - 1 AS start,
    end,
    Artifact_Type
  FROM 
    [_BLACKLISTED_TABLE_]) AS bl
ON
  seq.reference_name = bl.reference_name
WHERE 
  seq.start >= bl.start AND
  seq.end <= bl.end