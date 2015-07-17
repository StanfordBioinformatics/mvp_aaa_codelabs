# Publication Plots - QC

## Setup


Plot theme

```r
plot_theme = theme_minimal(base_size = 14, base_family = "Helvetica") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid = element_blank())
```

## QC Plots
Tables for sample level qc plots

```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.genome_calls_seq_qc",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.multi_sample_variants_seq_qc",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions",
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data")
```

#### Genotyping Concordance

```r
concordanceResult <- DisplayAndDispatchQuery("../sql/genotyping-concordance.sql",
                                             project=project,
                                             replacements=queryReplacements)
```

```
SELECT
  sample_id,
  calls_in_common,
  identical_calls,
  (identical_calls/calls_in_common) AS concordance
FROM (
  SELECT 
    sample_id,
    COUNT(seq_genotype) AS calls_in_common,
    SUM(IF(seq_genotype = gen_genotype, 1, 0)) AS identical_calls,
  FROM (
    SELECT
      seq.sample_id AS sample_id,
      seq.reference_name AS reference_name,
      seq.start AS start,
      seq.end AS end,
      seq.genotype AS seq_genotype,
      gen.genotype AS gen_genotype,
    FROM (
      SELECT
        sample_id,
        reference_name,
        start,
        end,
        genotype,
        bin,
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
          [va_aaa_pilot_data.genome_calls_seq_qc]
         #_WHERE_
        OMIT 
          call IF EVERY (call.genotype < 0)
        HAVING 
          num_alts <= 1
          AND reference_bases IN ('A','C','G','T')
          AND (alts IS null
            OR LENGTH(alts) <= 1)
        ),
        // Start javascript function
        // Input Columns
        call.call_set_name, reference_name, start, end, call.genotype,
        // Output Schema
        "[{name: 'sample_id', type: 'string'},
        {name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'end', type: 'integer'},
        {name: 'genotype', type: 'string'},
        {name: 'bin', type: 'integer'}]",
        // Function
        "function(r, emit) {
          for (c of r.call) {
            var binSize = 5000;
            var startBin = Math.floor(r.start/binSize);
            var endBin = Math.floor(r.end/binSize);
            var genotype = JSON.stringify(c.genotype.sort());
            for (var bin = startBin; bin <= endBin; bin++){
              emit({
                sample_id: c.call_set_name,
                reference_name: r.reference_name,
                start: r.start,
                end: r.end,
                genotype: genotype,
                bin: bin,
              })
            }
          }
        }")) AS seq
JOIN EACH (
  SELECT
    sample_id,
    reference_name,
    start,
    end,
    genotype,
    bin,
  FROM js(
    (SELECT
      call.call_set_name,
      reference_name,
      start,
      end,
      call.genotype,
    FROM
      [va_aaa_pilot_data.genotyping_data]
      OMIT call IF EVERY (call.genotype < 0)       
    ),
    // Start javascript function
    // Input Columns
    call.call_set_name, reference_name, start, end, call.genotype,
    // Output Schema
    "[{name: 'sample_id', type: 'string'},
    {name: 'reference_name', type: 'string'},
    {name: 'start', type: 'integer'},
    {name: 'end', type: 'integer'},
    {name: 'genotype', type: 'string'},
    {name: 'bin', type: 'integer'}]",
    // Function
    "function(r, emit) {
      for (c of r.call) {
        var binSize = 5000;
        var bin = Math.floor(r.start/binSize);
        var genotype = JSON.stringify(c.genotype.sort());
        var re = /\d+/;
        var chr = 'chr' + r.reference_name.match(re);
        emit({
          sample_id: c.call_set_name,
          reference_name: chr,
          start: r.start,
          end: r.end,
          genotype: genotype,
          bin: bin,
        })
      }
    }")) AS gen
ON
  seq.sample_id = gen.sample_id
  AND seq.reference_name = gen.reference_name
  AND seq.bin = gen.bin
WHERE
  seq.start <= gen.start
  AND seq.end >= gen.end )
GROUP BY 
  sample_id)
Running query:   RUNNING  2.5s
Running query:   RUNNING  3.1s
Running query:   RUNNING  3.8s
Running query:   RUNNING  4.4s
Running query:   RUNNING  5.0s
Running query:   RUNNING  5.6s
Running query:   RUNNING  6.2s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.5s
Running query:   RUNNING  8.1s
Running query:   RUNNING  8.8s
Running query:   RUNNING  9.4s
Running query:   RUNNING 10.0s
Running query:   RUNNING 10.6s
Running query:   RUNNING 11.2s
Running query:   RUNNING 11.8s
Running query:   RUNNING 12.4s
Running query:   RUNNING 13.0s
Running query:   RUNNING 13.6s
Running query:   RUNNING 14.2s
Running query:   RUNNING 14.9s
Running query:   RUNNING 15.5s
Running query:   RUNNING 16.2s
Running query:   RUNNING 16.8s
Running query:   RUNNING 17.5s
Running query:   RUNNING 18.1s
Running query:   RUNNING 18.7s
Running query:   RUNNING 19.3s
Running query:   RUNNING 19.9s
Running query:   RUNNING 20.5s
Running query:   RUNNING 21.1s
Running query:   RUNNING 21.7s
Running query:   RUNNING 22.3s
Running query:   RUNNING 22.9s
Running query:   RUNNING 23.5s
Running query:   RUNNING 24.1s
Running query:   RUNNING 24.8s
Running query:   RUNNING 25.4s
Running query:   RUNNING 26.0s
Running query:   RUNNING 26.6s
Running query:   RUNNING 27.2s
Running query:   RUNNING 27.9s
Running query:   RUNNING 28.5s
Running query:   RUNNING 29.1s
Running query:   RUNNING 29.7s
Running query:   RUNNING 30.3s
Running query:   RUNNING 30.9s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.2s
Running query:   RUNNING 32.8s
Running query:   RUNNING 33.4s
Running query:   RUNNING 34.0s
Running query:   RUNNING 34.6s
Running query:   RUNNING 35.2s
Running query:   RUNNING 35.9s
Running query:   RUNNING 36.5s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.8s
Running query:   RUNNING 38.4s
Running query:   RUNNING 39.0s
Running query:   RUNNING 39.6s
Running query:   RUNNING 40.2s
Running query:   RUNNING 40.8s
Running query:   RUNNING 41.5s
Running query:   RUNNING 42.1s
Running query:   RUNNING 42.7s
Running query:   RUNNING 43.3s
Running query:   RUNNING 44.0s
Running query:   RUNNING 44.6s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.8s
Running query:   RUNNING 46.4s
Running query:   RUNNING 47.0s
Running query:   RUNNING 47.6s
Running query:   RUNNING 48.2s
Running query:   RUNNING 48.9s
Running query:   RUNNING 49.5s
Running query:   RUNNING 50.1s
Running query:   RUNNING 50.7s
Running query:   RUNNING 51.3s
Running query:   RUNNING 51.9s
Running query:   RUNNING 52.5s
Running query:   RUNNING 53.1s
Running query:   RUNNING 53.8s
Running query:   RUNNING 54.4s
Running query:   RUNNING 55.0s
Running query:   RUNNING 55.6s
Running query:   RUNNING 56.2s
Running query:   RUNNING 56.8s
Running query:   RUNNING 57.5s
Running query:   RUNNING 58.1s
Running query:   RUNNING 58.8s
Running query:   RUNNING 59.4s
Running query:   RUNNING 60.0s
Running query:   RUNNING 60.6s
Running query:   RUNNING 61.2s
Running query:   RUNNING 61.9s
Running query:   RUNNING 62.5s
Running query:   RUNNING 63.1s
Running query:   RUNNING 63.8s
Running query:   RUNNING 64.4s
Running query:   RUNNING 65.0s
Running query:   RUNNING 65.6s
Running query:   RUNNING 66.2s
Running query:   RUNNING 66.8s
Running query:   RUNNING 67.5s
Running query:   RUNNING 68.1s
Running query:   RUNNING 68.7s
Running query:   RUNNING 69.3s
Running query:   RUNNING 69.9s
Running query:   RUNNING 70.5s
Running query:   RUNNING 71.1s
Running query:   RUNNING 71.7s
Running query:   RUNNING 72.3s
Running query:   RUNNING 72.9s
Running query:   RUNNING 73.5s
Running query:   RUNNING 74.1s
Running query:   RUNNING 74.7s
Running query:   RUNNING 75.3s
Running query:   RUNNING 76.0s
Running query:   RUNNING 76.6s
Running query:   RUNNING 77.2s
Running query:   RUNNING 77.8s
Running query:   RUNNING 78.4s
Running query:   RUNNING 79.0s
Running query:   RUNNING 79.6s
Running query:   RUNNING 80.2s
Running query:   RUNNING 80.9s
Running query:   RUNNING 81.5s
Running query:   RUNNING 82.1s
Running query:   RUNNING 82.8s
Running query:   RUNNING 83.4s
Running query:   RUNNING 84.2s
Running query:   RUNNING 84.9s
Running query:   RUNNING 85.5s
Running query:   RUNNING 86.1s
Running query:   RUNNING 86.7s
Running query:   RUNNING 87.3s
Running query:   RUNNING 88.0s
Running query:   RUNNING 88.6s
Running query:   RUNNING 89.2s
Running query:   RUNNING 89.9s
Running query:   RUNNING 90.5s
Running query:   RUNNING 91.1s
Running query:   RUNNING 91.7s
Running query:   RUNNING 92.4s
Running query:   RUNNING 93.0s
Running query:   RUNNING 93.6s
Running query:   RUNNING 94.3s
Running query:   RUNNING 95.0s
Running query:   RUNNING 95.6s
Running query:   RUNNING 96.2s
Running query:   RUNNING 97.1s
Running query:   RUNNING 97.9s
Running query:   RUNNING 98.5s
Running query:   RUNNING 99.1s
Running query:   RUNNING 99.7s
Running query:   RUNNING 100.3s
Running query:   RUNNING 100.9s
Running query:   RUNNING 101.6s
Running query:   RUNNING 102.2s
Running query:   RUNNING 102.8s
Running query:   RUNNING 103.5s
Running query:   RUNNING 104.1s
Running query:   RUNNING 104.7s
Running query:   RUNNING 105.3s
Running query:   RUNNING 105.9s
Running query:   RUNNING 106.6s
Running query:   RUNNING 107.2s
Running query:   RUNNING 107.9s
Running query:   RUNNING 108.5s
Running query:   RUNNING 109.1s
Running query:   RUNNING 109.7s
Running query:   RUNNING 110.3s
Running query:   RUNNING 111.0s
Running query:   RUNNING 111.6s
Running query:   RUNNING 112.2s
Running query:   RUNNING 112.8s
Running query:   RUNNING 113.4s
Running query:   RUNNING 114.1s
Running query:   RUNNING 114.7s
Running query:   RUNNING 115.3s
Running query:   RUNNING 115.9s
Running query:   RUNNING 116.5s
Running query:   RUNNING 117.1s
Running query:   RUNNING 117.7s
Running query:   RUNNING 118.3s
Running query:   RUNNING 118.9s
Running query:   RUNNING 119.5s
Running query:   RUNNING 120.1s
Running query:   RUNNING 120.8s
Running query:   RUNNING 121.4s
Running query:   RUNNING 122.0s
Running query:   RUNNING 122.6s
Running query:   RUNNING 123.2s
Running query:   RUNNING 123.8s
Running query:   RUNNING 124.4s
Running query:   RUNNING 125.0s
Running query:   RUNNING 125.7s
Running query:   RUNNING 126.3s
Running query:   RUNNING 126.9s
Running query:   RUNNING 127.5s
Running query:   RUNNING 128.2s
Running query:   RUNNING 128.8s
Running query:   RUNNING 129.4s
Running query:   RUNNING 130.0s
Running query:   RUNNING 130.6s
Running query:   RUNNING 131.2s
Running query:   RUNNING 131.8s
Running query:   RUNNING 132.4s
Running query:   RUNNING 133.0s
Running query:   RUNNING 133.6s
Running query:   RUNNING 134.3s
Running query:   RUNNING 134.9s
Running query:   RUNNING 135.5s
Running query:   RUNNING 136.1s
Running query:   RUNNING 136.7s
Running query:   RUNNING 137.3s
Running query:   RUNNING 138.0s
Running query:   RUNNING 138.6s
Running query:   RUNNING 139.2s
Running query:   RUNNING 139.8s
Running query:   RUNNING 140.4s
Running query:   RUNNING 141.1s
Running query:   RUNNING 141.7s
Running query:   RUNNING 142.3s
Running query:   RUNNING 143.0s
Running query:   RUNNING 143.6s
Running query:   RUNNING 144.2s
Running query:   RUNNING 144.8s
Running query:   RUNNING 145.4s
Running query:   RUNNING 146.0s
Running query:   RUNNING 146.6s
Running query:   RUNNING 147.2s
Running query:   RUNNING 147.8s
Running query:   RUNNING 148.5s
Running query:   RUNNING 149.1s
Running query:   RUNNING 149.7s
Running query:   RUNNING 150.3s
Running query:   RUNNING 150.9s
Running query:   RUNNING 151.5s
Running query:   RUNNING 152.2s
Running query:   RUNNING 152.8s
Running query:   RUNNING 153.4s
Running query:   RUNNING 154.0s
Running query:   RUNNING 154.6s
Running query:   RUNNING 155.3s
Running query:   RUNNING 155.9s
Running query:   RUNNING 156.5s
Running query:   RUNNING 157.1s
Running query:   RUNNING 157.7s
Running query:   RUNNING 158.3s
Running query:   RUNNING 158.9s
Running query:   RUNNING 159.5s
Running query:   RUNNING 160.2s
Running query:   RUNNING 160.9s
Running query:   RUNNING 161.5s
Running query:   RUNNING 162.1s
Running query:   RUNNING 162.8s
Running query:   RUNNING 163.4s
Running query:   RUNNING 164.0s
Running query:   RUNNING 164.6s
Running query:   RUNNING 165.2s
Running query:   RUNNING 165.9s
Running query:   RUNNING 166.5s
Running query:   RUNNING 167.1s
Running query:   RUNNING 167.7s
Running query:   RUNNING 168.4s
Running query:   RUNNING 169.0s
Running query:   RUNNING 169.6s
Running query:   RUNNING 170.2s
Running query:   RUNNING 170.8s
Running query:   RUNNING 171.4s
Running query:   RUNNING 172.0s
Running query:   RUNNING 172.6s
Running query:   RUNNING 173.3s
Running query:   RUNNING 173.9s
Running query:   RUNNING 174.5s
Running query:   RUNNING 175.1s
Running query:   RUNNING 175.8s
Running query:   RUNNING 176.4s
Running query:   RUNNING 177.0s
Running query:   RUNNING 177.6s
Running query:   RUNNING 178.2s
Running query:   RUNNING 178.9s
Running query:   RUNNING 179.5s
Running query:   RUNNING 180.1s
Running query:   RUNNING 180.8s
Running query:   RUNNING 181.4s
Running query:   RUNNING 182.1s
Running query:   RUNNING 182.7s
Running query:   RUNNING 183.3s
Running query:   RUNNING 183.9s
Running query:   RUNNING 184.6s
Running query:   RUNNING 185.2s
Running query:   RUNNING 185.8s
Running query:   RUNNING 186.4s
Running query:   RUNNING 187.0s
Running query:   RUNNING 187.6s
Running query:   RUNNING 188.3s
Running query:   RUNNING 188.9s
Running query:   RUNNING 189.5s
Running query:   RUNNING 190.1s
Running query:   RUNNING 190.7s
Running query:   RUNNING 191.3s
Running query:   RUNNING 191.9s
Running query:   RUNNING 192.6s
Running query:   RUNNING 193.2s
Running query:   RUNNING 193.8s
Running query:   RUNNING 194.4s
Running query:   RUNNING 195.0s
Running query:   RUNNING 195.6s
Running query:   RUNNING 196.2s
Running query:   RUNNING 196.9s
Running query:   RUNNING 197.5s
Running query:   RUNNING 198.1s
Running query:   RUNNING 198.8s
Running query:   RUNNING 199.4s
Running query:   RUNNING 200.0s
Running query:   RUNNING 200.6s
Running query:   RUNNING 201.2s
Running query:   RUNNING 201.8s
Running query:   RUNNING 202.4s
Running query:   RUNNING 203.0s
Running query:   RUNNING 203.7s
Running query:   RUNNING 204.3s
Running query:   RUNNING 205.0s
Running query:   RUNNING 205.6s
Running query:   RUNNING 206.2s
Running query:   RUNNING 206.8s
Running query:   RUNNING 207.4s
Running query:   RUNNING 208.0s
Running query:   RUNNING 208.7s
Running query:   RUNNING 209.3s
Running query:   RUNNING 209.9s
Running query:   RUNNING 210.5s
Running query:   RUNNING 211.1s
Running query:   RUNNING 211.7s
Running query:   RUNNING 212.3s
Running query:   RUNNING 213.0s
Running query:   RUNNING 213.6s
Running query:   RUNNING 214.2s
Running query:   RUNNING 214.8s
Running query:   RUNNING 215.4s
Running query:   RUNNING 216.1s
Running query:   RUNNING 216.7s
Running query:   RUNNING 217.3s
Running query:   RUNNING 217.9s
Running query:   RUNNING 218.5s
Running query:   RUNNING 219.1s
Running query:   RUNNING 219.8s
Running query:   RUNNING 220.4s
Running query:   RUNNING 221.0s
Running query:   RUNNING 221.6s
Running query:   RUNNING 222.2s
Running query:   RUNNING 222.9s
Running query:   RUNNING 223.5s
Running query:   RUNNING 224.1s
Running query:   RUNNING 224.7s
Running query:   RUNNING 225.4s
Running query:   RUNNING 226.2s
Running query:   RUNNING 226.8s
Running query:   RUNNING 227.5s
Running query:   RUNNING 228.2s
Running query:   RUNNING 228.8s
Running query:   RUNNING 229.4s
Running query:   RUNNING 230.0s
Running query:   RUNNING 230.6s
Running query:   RUNNING 231.2s
Running query:   RUNNING 231.8s
Running query:   RUNNING 232.4s
Running query:   RUNNING 233.1s
Running query:   RUNNING 233.7s
Running query:   RUNNING 234.3s
Running query:   RUNNING 234.9s
Running query:   RUNNING 235.5s
Running query:   RUNNING 236.1s
Running query:   RUNNING 236.7s
Running query:   RUNNING 237.3s
Running query:   RUNNING 238.0s
Running query:   RUNNING 238.6s
Running query:   RUNNING 239.2s
Running query:   RUNNING 239.8s
Running query:   RUNNING 240.4s
Running query:   RUNNING 241.0s
Running query:   RUNNING 241.7s
Running query:   RUNNING 242.3s
Running query:   RUNNING 242.9s
Running query:   RUNNING 243.5s
Running query:   RUNNING 244.1s
Running query:   RUNNING 244.7s
Running query:   RUNNING 245.3s
Running query:   RUNNING 246.0s
Running query:   RUNNING 246.6s
Running query:   RUNNING 247.2s
Running query:   RUNNING 247.9s
Running query:   RUNNING 248.6s
Running query:   RUNNING 249.2s
Running query:   RUNNING 249.9s
Running query:   RUNNING 250.5s
Running query:   RUNNING 251.1s
Running query:   RUNNING 251.7s
Running query:   RUNNING 252.5s
Running query:   RUNNING 253.1s
Running query:   RUNNING 253.7s
Running query:   RUNNING 254.3s
Running query:   RUNNING 254.9s
Running query:   RUNNING 255.5s
Running query:   RUNNING 256.1s
Running query:   RUNNING 256.7s
Running query:   RUNNING 257.4s
Running query:   RUNNING 258.0s
Running query:   RUNNING 258.6s
Running query:   RUNNING 259.2s
Running query:   RUNNING 259.8s
Running query:   RUNNING 260.4s
Running query:   RUNNING 261.0s
Running query:   RUNNING 261.7s
Running query:   RUNNING 262.3s
Running query:   RUNNING 263.0s
Running query:   RUNNING 263.6s
Running query:   RUNNING 264.2s
Running query:   RUNNING 265.0s
Running query:   RUNNING 265.6s
Running query:   RUNNING 266.2s
Running query:   RUNNING 266.9s
Running query:   RUNNING 267.5s
Running query:   RUNNING 268.1s
Running query:   RUNNING 268.7s
Running query:   RUNNING 269.3s
Running query:   RUNNING 269.9s
Running query:   RUNNING 270.5s
Running query:   RUNNING 271.1s
Running query:   RUNNING 271.8s
Running query:   RUNNING 272.4s
Running query:   RUNNING 273.0s
Running query:   RUNNING 273.6s
Running query:   RUNNING 274.2s
Running query:   RUNNING 274.9s
Running query:   RUNNING 275.5s
Running query:   RUNNING 276.1s
Running query:   RUNNING 276.8s
Running query:   RUNNING 277.4s
Running query:   RUNNING 278.0s
Running query:   RUNNING 278.6s
Running query:   RUNNING 279.2s
Running query:   RUNNING 279.8s
Running query:   RUNNING 280.5s
Running query:   RUNNING 281.1s
Running query:   RUNNING 281.8s
Running query:   RUNNING 282.4s
Running query:   RUNNING 283.1s
Running query:   RUNNING 283.7s
Running query:   RUNNING 284.3s
Running query:   RUNNING 284.9s
Running query:   RUNNING 285.5s
Running query:   RUNNING 286.2s
Running query:   RUNNING 286.9s
Running query:   RUNNING 287.6s
Running query:   RUNNING 288.2s
Running query:   RUNNING 288.8s
Running query:   RUNNING 289.5s
Running query:   RUNNING 290.1s
Running query:   RUNNING 290.7s
Running query:   RUNNING 291.3s
Running query:   RUNNING 291.9s
Running query:   RUNNING 292.5s
Running query:   RUNNING 293.1s
Running query:   RUNNING 293.8s
Running query:   RUNNING 294.4s
Running query:   RUNNING 295.0s
Running query:   RUNNING 295.6s
Running query:   RUNNING 296.2s
Running query:   RUNNING 296.9s
Running query:   RUNNING 297.5s
Running query:   RUNNING 298.1s
Running query:   RUNNING 298.7s
Running query:   RUNNING 299.3s
Running query:   RUNNING 299.9s
Running query:   RUNNING 300.5s
Running query:   RUNNING 301.1s
Running query:   RUNNING 301.7s
Running query:   RUNNING 302.3s
Running query:   RUNNING 302.9s
Running query:   RUNNING 303.5s
Running query:   RUNNING 304.2s
Running query:   RUNNING 304.8s
Running query:   RUNNING 305.4s
Running query:   RUNNING 306.0s
Running query:   RUNNING 306.6s
Running query:   RUNNING 307.3s
Running query:   RUNNING 307.9s
Running query:   RUNNING 308.5s
Running query:   RUNNING 309.1s
Running query:   RUNNING 309.8s
Running query:   RUNNING 310.4s
Running query:   RUNNING 311.0s
Running query:   RUNNING 311.6s
Running query:   RUNNING 312.2s
Running query:   RUNNING 312.8s
Running query:   RUNNING 313.4s
Running query:   RUNNING 314.0s
Running query:   RUNNING 314.6s
Running query:   RUNNING 315.3s
Running query:   RUNNING 315.9s
Running query:   RUNNING 316.5s
Running query:   RUNNING 317.1s
Running query:   RUNNING 317.7s
Running query:   RUNNING 318.4s
Running query:   RUNNING 319.0s
Running query:   RUNNING 319.6s
Running query:   RUNNING 320.2s
Running query:   RUNNING 320.8s
Running query:   RUNNING 321.4s
Running query:   RUNNING 322.0s
Running query:   RUNNING 322.6s
Running query:   RUNNING 323.3s
Running query:   RUNNING 323.9s
Running query:   RUNNING 324.5s
Running query:   RUNNING 325.1s
Running query:   RUNNING 325.7s
Running query:   RUNNING 326.3s
Running query:   RUNNING 326.9s
Running query:   RUNNING 327.5s
Running query:   RUNNING 328.1s
Running query:   RUNNING 328.7s
Running query:   RUNNING 329.4s
Running query:   RUNNING 330.0s
Running query:   RUNNING 330.6s
Running query:   RUNNING 331.2s
Running query:   RUNNING 331.8s
Running query:   RUNNING 332.4s
Running query:   RUNNING 333.0s
Running query:   RUNNING 333.7s
Running query:   RUNNING 334.3s
Running query:   RUNNING 334.9s
Running query:   RUNNING 335.6s
Running query:   RUNNING 336.2s
Running query:   RUNNING 336.8s
Running query:   RUNNING 337.4s
Running query:   RUNNING 338.0s
Running query:   RUNNING 338.6s
Running query:   RUNNING 339.2s
Running query:   RUNNING 339.8s
Running query:   RUNNING 340.5s
Running query:   RUNNING 341.1s
Running query:   RUNNING 341.7s
Running query:   RUNNING 342.3s
Running query:   RUNNING 342.9s
Running query:   RUNNING 343.5s
Running query:   RUNNING 344.1s
Running query:   RUNNING 344.8s
Running query:   RUNNING 345.4s
Running query:   RUNNING 346.0s
Running query:   RUNNING 346.6s
Running query:   RUNNING 347.2s
Running query:   RUNNING 347.9s
Running query:   RUNNING 348.5s
Running query:   RUNNING 349.1s
Running query:   RUNNING 349.7s
Running query:   RUNNING 350.3s
Running query:   RUNNING 350.9s
Running query:   RUNNING 351.5s
Running query:   RUNNING 352.1s
Running query:   RUNNING 352.7s
Running query:   RUNNING 353.4s
Running query:   RUNNING 354.0s
Running query:   RUNNING 354.6s
Running query:   RUNNING 355.2s
Running query:   RUNNING 355.8s
Running query:   RUNNING 356.4s
Running query:   RUNNING 357.0s
Running query:   RUNNING 357.6s
Running query:   RUNNING 358.2s
Running query:   RUNNING 358.9s
Running query:   RUNNING 359.5s
Running query:   RUNNING 360.1s
```

```r
plate = substr(concordanceResult$sample_id, 1, 9)
concordanceResult = cbind(concordanceResult, plate)
```


```r
ggplot(concordanceResult) +
  geom_point(aes(x=sample_id, y=concordance, color=plate), size=4, alpha=0.5) +
  xlab("Sample") +
  ylab("Concordance") +
  ggtitle("Concordance with Genotyping Data") +
  scale_colour_brewer(name="Library Prep Plate", palette="Set1") +
  plot_theme + 
  theme(axis.text.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.2, 0.8))
```

<img src="figure/concordance-publication-1.png" title="plot of chunk concordance-publication" alt="plot of chunk concordance-publication" style="display: block; margin: auto;" />

####Batch Effect PCA

```r
pcaFile = '/Users/gmcinnes/data/pca-all-genomes-all-references-no-tissue.tsv'
pcaResult = read.table(pcaFile)
names(pcaResult) = c('sample_id','pc1', 'pc2', 'something')
plate = substr(pcaResult$sample_id, 1, 9)
pcaResult = cbind(pcaResult, plate)
```


```r
ggplot(pcaResult, aes(pc1, pc2, color=plate)) + 
  geom_point(size=4, alpha=0.5) +
  ggtitle("Principal Component Analysis") +
  xlab("Principal Component 1") + 
  ylab("Principal Component 2") +
  scale_colour_brewer(name="Library Prep Plate", palette="Set1") +
  plot_theme +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position = c(0.85,0.75))
```

<img src="figure/pca-batch-effect-publication-1.png" title="plot of chunk pca-batch-effect-publication" alt="plot of chunk pca-batch-effect-publication" style="display: block; margin: auto;" />

#### IBS

```r
ibsDataFlowFilename = '../data/ibs-460.tsv'
ReadIBSFile <- function(ibsFilename, header=FALSE, rowNames=NULL) {
  ibsData <- read.table(file=ibsFilename, header=header,
                        row.names=rowNames, stringsAsFactors=FALSE)
  return (ibsData)
}
ibsDataflowData <- ReadIBSFile(ibsDataFlowFilename)

ColumnNames <- function(ibsData) { 
  if(3 == ncol(ibsData)) {
    colnames(ibsData) <- c("sample1", "sample2", "ibsScore")
  } else {
    colnames(ibsData) <- c("sample1", "sample2", "ibsScore", "similar", "observed")
  }
}
colnames(ibsDataflowData) <- ColumnNames(ibsDataflowData)

MakeIBSDataSymmetric <- function(ibsData) {
  ibsPairsMirrored <- data.frame(sample1=ibsData$sample2,
                                 sample2=ibsData$sample1,
                                 ibsScore=ibsData$ibsScore)
  ibsData <- rbind(ibsData[,1:3], ibsPairsMirrored)
}
ibsDataflowData <- MakeIBSDataSymmetric(ibsDataflowData)

ExcludeDiagonal <- function(ibsData) {
  ibsData <- filter(ibsData, ibsData$sample1 != ibsData$sample2)
  return (ibsData)
}
ibsDataflowDataSample <- ExcludeDiagonal(ibsDataflowData)


SampleIBSMatrix <- function(ibsData, sampleSize=50) {
  individuals <- unique(ibsData$sample1)
  sample <- sample(individuals, sampleSize)
  # In order to highlight related samples we make sure these are included in the subset
  related = c('LP6005051-DNA_D09',
              'LP6005051-DNA_E02', 
              'LP6005692-DNA_D05', 
              'LP6005243-DNA_E10',
              'LP6005144-DNA_A02',
              "LP6005051-DNA_D04",
              "LP6005243-DNA_H03",
              "LP6005144-DNA_D04",
              "LP6005692-DNA_E10",
              "LP6005692-DNA_G09")
  sample = c(sample, related)
  ibsData <- subset(ibsData, ibsData$sample1 %in% sample)
  ibsData <- subset(ibsData, ibsData$sample2 %in% sample)
  return (ibsData)
}
ibsDataflowDataSubset <- SampleIBSMatrix(ibsDataflowDataSample)

DrawHeatMap <- function(ibsData) {
  p <- ggplot(data=ibsData, aes(x=sample1, y=sample2)) +
    plot_theme +
    theme(axis.ticks=element_blank(), 
          axis.text=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_tile(aes(fill=ibsScore), colour="white") +
    scale_fill_gradient(low="white", high="steelblue", na.value="black",
                        guide=guide_colourbar(title= "IBS Score")) +
    labs(list(title="Identity By State (IBS) Heat Map",
              x="Sample", y="Sample")) 
  p 
}
```


```r
DrawHeatMap(ibsDataflowDataSubset)
```

<img src="figure/ibs-publication-1.png" title="plot of chunk ibs-publication" alt="plot of chunk ibs-publication" style="display: block; margin: auto;" />


#### Sex Inference

```r
result <- DisplayAndDispatchQuery("../sql/gender-check.sql",
                                  project=project,
                                  replacements=queryReplacements)
```

```
# Compute the the homozygous and heterozygous variant counts for each individual
# within chromosome X to help determine whether the gender phenotype value is
# correct for each individual.
SELECT
  call.call_set_name,
  ROUND((het_RA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_het_alt_in_snvs,
  ROUND((hom_AA_count/(hom_AA_count + het_RA_count))*1000)/1000 AS perct_hom_alt_in_snvs,
  (hom_AA_count + het_RA_count + hom_RR_count) AS all_callable_sites,
  hom_AA_count,
  het_RA_count,
  hom_RR_count,
  (hom_AA_count + het_RA_count) AS all_snvs,
FROM
  (
  SELECT
    call.call_set_name,
    SUM(0 = first_allele
      AND 0 = second_allele) AS hom_RR_count,
    SUM(first_allele = second_allele AND first_allele > 0) AS hom_AA_count,
    SUM((first_allele != second_allele OR second_allele IS NULL)
      AND (first_allele > 0 OR second_allele > 0)) AS het_RA_count
  FROM (
    SELECT
      reference_bases,
      GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.call_set_name,
      NTH(1, call.genotype) WITHIN call AS first_allele,
      NTH(2, call.genotype) WITHIN call AS second_allele,
    FROM
      [va_aaa_pilot_data.multi_sample_variants_seq_qc]
    WHERE
      reference_name = 'chrX'
      AND start NOT BETWEEN 59999 AND 2699519
      AND start NOT BETWEEN 154931042 AND 155260559
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T')
      )
  GROUP BY
    call.call_set_name)
ORDER BY
  call.call_set_name
```

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender), size=4, alpha=0.5) +
  xlab("Sample") +
  ylab("Heterozygosity Rate") +
  ggtitle("Heterozygosity Rate on the X Chromosome") +
  scale_colour_brewer(palette="Set1", name="Gender") +
  plot_theme +
  theme(axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = c(0.85, 0.55))
```

<img src="figure/sex-inference-publication-1.png" title="plot of chunk sex-inference-publication" alt="plot of chunk sex-inference-publication" style="display: block; margin: auto;" />


Tables for sample level qc plots

```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.genome_calls_seq_qc",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.multi_sample_variants_sample_qc",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions",
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data")
```

#### Ti/Tv By Genomic Window

```r
titvWindowResults <- DisplayAndDispatchQuery("../sql/ti-tv-ratio.sql",
                                             project=project,
                                             replacements=c("#_WHERE_"="WHERE reference_name = 'chr1'",
                                                            "_WINDOW_SIZE_"="100000",
                                                            queryReplacements))
```

```
# Compute the Ti/Tv ratio for variants within genomic region windows.
SELECT
  reference_name,
  window * 100000 AS window_start,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_window,
FROM (
  SELECT
    reference_name,
    window,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    COUNT(mutation) AS num_variants_in_window
  FROM (
    SELECT
      reference_name,
      INTEGER(FLOOR(start / 100000)) AS window,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    FROM
      [va_aaa_pilot_data.genome_calls_seq_qc]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    WHERE reference_name = 'chr1'
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
  GROUP BY
    reference_name,
    window)
ORDER BY
  window_start
```


```r
ggplot(titvWindowResults, aes(x=window_start, y=titv)) +
  geom_point() +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  scale_x_continuous(expand = c(0, 0), labels=comma) +
  ggtitle("Ti/Tv by 100,000 base pair\nwindows on Chromosome 1") +
  plot_theme
```

<img src="figure/titv-by-genomic-window-publication-1.png" title="plot of chunk titv-by-genomic-window-publication" alt="plot of chunk titv-by-genomic-window-publication" style="display: block; margin: auto;" />

#### Ti/Tv By Depth

```r
query <- "../sql/ti-tv-by-depth.sql"
titv <- DisplayAndDispatchQuery(query,
                                project=project,
                                replacements=c(queryReplacements))
```

```
SELECT
  call.call_set_name,
  (transitions/transversions) AS titv_ratio,
  average_depth,
FROM (
  SELECT
    call.call_set_name,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    ROUND(AVG(call.DP)) AS average_depth,
  FROM (

    SELECT
      call.call_set_name,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
      call.DP
    FROM (
      SELECT
        call.call_set_name,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        call.genotype,
        call.DP,
      FROM
        [va_aaa_pilot_data.multi_sample_variants_sample_qc]
      # Optionally add clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_  
      )
    WHERE
      call.DP is not null
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
    GROUP BY 
      call.call_set_name,
      call.DP,)
WHERE
  transversions > 0
GROUP BY
  call.call_set_name,
  titv_ratio,
  average_depth,

Retrieving data:  3.1s
Retrieving data:  4.3s
Retrieving data:  5.6s
Retrieving data:  6.9s
Retrieving data:  8.2s
Retrieving data:  9.5s
Retrieving data: 10.8s
Retrieving data: 11.9s
Retrieving data: 13.2s
```


```r
ggplot(titv, aes(x=average_depth, y=titv_ratio, color=call_call_set_name)) + 
  geom_point(size=2) +
  ggtitle("Ti/Tv By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  plot_theme + 
  theme(legend.position="none")
```

<img src="figure/titv-by-depth-publication-1.png" title="plot of chunk titv-by-depth-publication" alt="plot of chunk titv-by-depth-publication" style="display: block; margin: auto;" />


