#SELECT
#  COUNT(start) AS REF_COUNT
#FROM (
SELECT
  reference_name,
  start,
  end,
    FROM js(
      (SELECT
reference_name,
bin,
GROUP_CONCAT(joined) AS blocks
FROM (
SELECT
  reference_name,
  start,
  end,
  bin,
  CONCAT(STRING(start), ":", STRING(end)) AS joined
FROM (
SELECT
  call.call_set_name,
  reference_name,
  start,
  end,
  FLOOR(start/5000) AS bin,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alts,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
FROM
  va_aaa_pilot_data.genome_calls_seq_qc
OMIT RECORD IF EVERY(alternate_bases IS NOT NULL)
  OR SOME(call.genotype != 0)
LIMIT 1000)
ORDER BY
  reference_name,
  start,
  end,)
GROUP BY
  reference_name,
  bin),
      reference_name, blocks, bin,
      # Find overlaps
      "[{name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'end', type: 'integer'}]",
       "function(r, emit) {
            //var blocks = r.blocks;
            var ranges = r.blocks.split(",");
            var last_start = 0;
            var last_end = 0;
            for	(var index = 0; index < ranges.length; index++) {
              var positions = ranges[index].split(":");
              var start = positions[0];
              var end = positions[1];
              if (start > last_end) {
                emit({
                  reference_name: r.reference_name,
                  start: last_start,
                  end: last_end,
                })
                last_start = start;
                last_end = end;
              }
              else {
                last_end = end;
              }
            }
        }")
#GROUP EACH BY
#  reference_name,
#  start,
#  end)
