# Get all the blacklisted variants in the dataset
SELECT
  seq.reference_name,
  seq.start,
  seq.end,
  bl.Artifact_Type
FROM (
  SELECT
    reference_name,
    start,
    end,
  FROM
    [_THE_TABLE_]) as seq
JOIN (
  SELECT 
    reference_name,
    start - 1 AS start,
    end,
    Artifact_Type
  FROM 
    [resources.blacklisted_positions]) AS bl
ON
  seq.reference_name = bl.reference_name
WHERE 
  seq.start >= bl.start AND
  seq.end <= bl.end