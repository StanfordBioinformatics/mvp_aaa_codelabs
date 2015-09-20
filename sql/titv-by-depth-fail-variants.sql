# Select variant ids for samples at a specified depth

SELECT
  variant_id,
  call.call_set_name AS sample_id,
  "titv_by_depth" AS failure_reason,
FROM
  [_THE_EXPANDED_TABLE_]
WHERE
  # A list of samples and the depths of interest
  # Example:
  #   (call.call_set_name = 'LP6005692-DNA_B08'
  #      AND call.DP = 112)
  #    OR (call.call_set_name = 'LP6005692-DNA_C10'
  #      AND call.DP = 125)
  _SAMPLE_DEPTH_
#_ORDER_BY_
