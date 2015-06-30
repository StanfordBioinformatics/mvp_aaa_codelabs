# Publication Plots - QC

## Setup



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
Running query:   RUNNING  2.8s
Running query:   RUNNING  3.6s
Running query:   RUNNING  5.1s
Running query:   RUNNING  5.8s
Running query:   RUNNING  6.4s
Running query:   RUNNING  7.2s
Running query:   RUNNING  7.9s
Running query:   RUNNING  8.6s
Running query:   RUNNING  9.3s
Running query:   RUNNING 10.0s
Running query:   RUNNING 10.7s
Running query:   RUNNING 11.4s
Running query:   RUNNING 12.5s
Running query:   RUNNING 13.2s
Running query:   RUNNING 14.0s
Running query:   RUNNING 14.7s
Running query:   RUNNING 15.4s
Running query:   RUNNING 16.2s
Running query:   RUNNING 17.2s
Running query:   RUNNING 18.7s
Running query:   RUNNING 19.5s
Running query:   RUNNING 20.2s
Running query:   RUNNING 21.7s
Running query:   RUNNING 22.6s
Running query:   RUNNING 23.3s
Running query:   RUNNING 24.2s
Running query:   RUNNING 24.9s
Running query:   RUNNING 26.0s
Running query:   RUNNING 26.8s
Running query:   RUNNING 27.6s
Running query:   RUNNING 29.2s
Running query:   RUNNING 30.0s
Running query:   RUNNING 30.7s
Running query:   RUNNING 31.5s
Running query:   RUNNING 32.2s
Running query:   RUNNING 33.8s
Running query:   RUNNING 37.0s
Running query:   RUNNING 39.3s
Running query:   RUNNING 40.9s
Running query:   RUNNING 41.7s
Running query:   RUNNING 42.4s
Running query:   RUNNING 43.1s
Running query:   RUNNING 43.8s
Running query:   RUNNING 44.5s
Running query:   RUNNING 45.2s
Running query:   RUNNING 45.9s
Running query:   RUNNING 46.7s
Running query:   RUNNING 47.4s
Running query:   RUNNING 48.1s
Running query:   RUNNING 48.8s
Running query:   RUNNING 49.6s
Running query:   RUNNING 50.3s
Running query:   RUNNING 51.0s
Running query:   RUNNING 51.7s
Running query:   RUNNING 52.5s
Running query:   RUNNING 53.3s
Running query:   RUNNING 54.0s
Running query:   RUNNING 54.7s
Running query:   RUNNING 55.5s
Running query:   RUNNING 56.2s
Running query:   RUNNING 56.9s
Running query:   RUNNING 57.6s
Running query:   RUNNING 58.4s
Running query:   RUNNING 59.1s
Running query:   RUNNING 59.8s
Running query:   RUNNING 60.6s
Running query:   RUNNING 61.3s
Running query:   RUNNING 62.6s
Running query:   RUNNING 63.5s
Running query:   RUNNING 67.3s
Running query:   RUNNING 68.4s
Running query:   RUNNING 69.1s
Running query:   RUNNING 70.0s
Running query:   RUNNING 71.8s
Running query:   RUNNING 72.9s
Running query:   RUNNING 74.3s
Running query:   RUNNING 75.2s
Running query:   RUNNING 76.7s
Running query:   RUNNING 77.6s
Running query:   RUNNING 78.2s
Running query:   RUNNING 78.9s
Running query:   RUNNING 79.6s
Running query:   RUNNING 80.3s
Running query:   RUNNING 81.1s
Running query:   RUNNING 81.8s
Running query:   RUNNING 82.5s
Running query:   RUNNING 83.2s
Running query:   RUNNING 83.9s
Running query:   RUNNING 84.7s
Running query:   RUNNING 85.4s
Running query:   RUNNING 86.1s
Running query:   RUNNING 86.8s
Running query:   RUNNING 87.6s
Running query:   RUNNING 88.3s
Running query:   RUNNING 89.0s
Running query:   RUNNING 89.7s
Running query:   RUNNING 90.4s
Running query:   RUNNING 91.1s
Running query:   RUNNING 91.9s
Running query:   RUNNING 92.6s
Running query:   RUNNING 93.3s
Running query:   RUNNING 94.0s
Running query:   RUNNING 94.7s
Running query:   RUNNING 95.4s
Running query:   RUNNING 96.2s
Running query:   RUNNING 96.9s
Running query:   RUNNING 97.6s
Running query:   RUNNING 98.3s
Running query:   RUNNING 99.0s
Running query:   RUNNING 99.8s
Running query:   RUNNING 102.0s
Running query:   RUNNING 103.5s
Running query:   RUNNING 104.7s
Running query:   RUNNING 105.5s
Running query:   RUNNING 106.3s
Running query:   RUNNING 107.1s
Running query:   RUNNING 108.2s
Running query:   RUNNING 109.0s
Running query:   RUNNING 109.9s
Running query:   RUNNING 110.6s
Running query:   RUNNING 111.3s
Running query:   RUNNING 112.2s
Running query:   RUNNING 112.9s
Running query:   RUNNING 114.1s
Running query:   RUNNING 114.7s
Running query:   RUNNING 115.5s
Running query:   RUNNING 116.2s
Running query:   RUNNING 116.9s
Running query:   RUNNING 117.6s
Running query:   RUNNING 118.3s
Running query:   RUNNING 119.0s
Running query:   RUNNING 119.8s
Running query:   RUNNING 120.7s
Running query:   RUNNING 121.5s
Running query:   RUNNING 122.2s
Running query:   RUNNING 122.9s
Running query:   RUNNING 123.6s
Running query:   RUNNING 124.4s
Running query:   RUNNING 125.1s
Running query:   RUNNING 125.8s
Running query:   RUNNING 126.5s
Running query:   RUNNING 127.2s
Running query:   RUNNING 128.0s
Running query:   RUNNING 128.7s
Running query:   RUNNING 129.4s
Running query:   RUNNING 130.2s
Running query:   RUNNING 130.9s
Running query:   RUNNING 131.6s
Running query:   RUNNING 132.3s
Running query:   RUNNING 133.4s
Running query:   RUNNING 134.2s
Running query:   RUNNING 134.9s
Running query:   RUNNING 135.6s
Running query:   RUNNING 136.4s
Running query:   RUNNING 137.1s
Running query:   RUNNING 137.9s
Running query:   RUNNING 138.6s
Running query:   RUNNING 139.3s
Running query:   RUNNING 140.1s
Running query:   RUNNING 140.8s
Running query:   RUNNING 141.5s
Running query:   RUNNING 142.2s
Running query:   RUNNING 143.0s
Running query:   RUNNING 143.7s
Running query:   RUNNING 144.4s
Running query:   RUNNING 145.5s
Running query:   RUNNING 146.2s
Running query:   RUNNING 147.0s
Running query:   RUNNING 147.7s
Running query:   RUNNING 148.5s
Running query:   RUNNING 149.2s
Running query:   RUNNING 149.9s
Running query:   RUNNING 150.8s
Running query:   RUNNING 151.6s
Running query:   RUNNING 152.7s
Running query:   RUNNING 153.4s
Running query:   RUNNING 154.1s
Running query:   RUNNING 154.8s
Running query:   RUNNING 155.5s
Running query:   RUNNING 156.2s
Running query:   RUNNING 157.1s
Running query:   RUNNING 157.8s
Running query:   RUNNING 158.6s
Running query:   RUNNING 159.4s
Running query:   RUNNING 160.5s
Running query:   RUNNING 161.3s
Running query:   RUNNING 162.0s
Running query:   RUNNING 162.8s
Running query:   RUNNING 163.5s
Running query:   RUNNING 164.7s
Running query:   RUNNING 165.5s
Running query:   RUNNING 166.3s
Running query:   RUNNING 167.6s
Running query:   RUNNING 168.3s
Running query:   RUNNING 169.0s
Running query:   RUNNING 170.3s
Running query:   RUNNING 170.9s
Running query:   RUNNING 171.6s
Running query:   RUNNING 172.3s
Running query:   RUNNING 173.0s
Running query:   RUNNING 173.8s
Running query:   RUNNING 174.5s
Running query:   RUNNING 175.2s
Running query:   RUNNING 175.9s
Running query:   RUNNING 176.6s
Running query:   RUNNING 177.4s
Running query:   RUNNING 178.1s
Running query:   RUNNING 178.8s
Running query:   RUNNING 179.5s
Running query:   RUNNING 180.2s
Running query:   RUNNING 180.9s
Running query:   RUNNING 181.7s
Running query:   RUNNING 182.8s
Running query:   RUNNING 183.5s
Running query:   RUNNING 184.2s
Running query:   RUNNING 184.9s
Running query:   RUNNING 185.7s
Running query:   RUNNING 186.4s
Running query:   RUNNING 187.1s
Running query:   RUNNING 187.8s
Running query:   RUNNING 188.5s
Running query:   RUNNING 189.2s
Running query:   RUNNING 190.0s
Running query:   RUNNING 190.8s
Running query:   RUNNING 191.5s
Running query:   RUNNING 192.2s
Running query:   RUNNING 192.9s
Running query:   RUNNING 194.1s
Running query:   RUNNING 194.9s
Running query:   RUNNING 195.6s
Running query:   RUNNING 197.7s
Running query:   RUNNING 198.5s
Running query:   RUNNING 199.2s
Running query:   RUNNING 200.0s
Running query:   RUNNING 201.4s
Running query:   RUNNING 202.7s
Running query:   RUNNING 203.5s
Running query:   RUNNING 204.2s
Running query:   RUNNING 204.9s
Running query:   RUNNING 205.7s
Running query:   RUNNING 206.4s
Running query:   RUNNING 207.2s
Running query:   RUNNING 208.2s
Running query:   RUNNING 208.9s
Running query:   RUNNING 209.6s
Running query:   RUNNING 210.5s
Running query:   RUNNING 211.2s
Running query:   RUNNING 211.9s
Running query:   RUNNING 212.6s
Running query:   RUNNING 213.4s
Running query:   RUNNING 214.1s
Running query:   RUNNING 214.8s
Running query:   RUNNING 215.6s
Running query:   RUNNING 216.3s
Running query:   RUNNING 217.0s
Running query:   RUNNING 217.7s
Running query:   RUNNING 218.5s
Running query:   RUNNING 219.2s
Running query:   RUNNING 219.9s
Running query:   RUNNING 220.8s
Running query:   RUNNING 221.5s
Running query:   RUNNING 222.3s
Running query:   RUNNING 223.0s
Running query:   RUNNING 223.8s
Running query:   RUNNING 224.5s
Running query:   RUNNING 225.2s
Running query:   RUNNING 225.9s
Running query:   RUNNING 226.6s
Running query:   RUNNING 227.4s
Running query:   RUNNING 228.1s
Running query:   RUNNING 228.8s
Running query:   RUNNING 229.5s
Running query:   RUNNING 230.6s
Running query:   RUNNING 231.3s
Running query:   RUNNING 232.1s
Running query:   RUNNING 232.8s
Running query:   RUNNING 234.7s
Running query:   RUNNING 235.4s
Running query:   RUNNING 236.2s
Running query:   RUNNING 237.8s
Running query:   RUNNING 238.5s
Running query:   RUNNING 239.3s
Running query:   RUNNING 240.2s
Running query:   RUNNING 241.4s
Running query:   RUNNING 242.0s
Running query:   RUNNING 242.8s
Running query:   RUNNING 243.5s
Running query:   RUNNING 244.2s
Running query:   RUNNING 244.9s
Running query:   RUNNING 245.6s
Running query:   RUNNING 246.4s
Running query:   RUNNING 247.1s
Running query:   RUNNING 247.8s
Running query:   RUNNING 248.5s
Running query:   RUNNING 249.2s
Running query:   RUNNING 249.9s
Running query:   RUNNING 250.6s
Running query:   RUNNING 251.3s
Running query:   RUNNING 252.1s
Running query:   RUNNING 252.8s
Running query:   RUNNING 253.5s
Running query:   RUNNING 254.2s
Running query:   RUNNING 255.0s
Running query:   RUNNING 255.9s
Running query:   RUNNING 256.6s
Running query:   RUNNING 257.3s
Running query:   RUNNING 258.4s
Running query:   RUNNING 259.1s
Running query:   RUNNING 259.9s
Running query:   RUNNING 260.6s
Running query:   RUNNING 261.3s
Running query:   RUNNING 262.0s
Running query:   RUNNING 262.7s
Running query:   RUNNING 263.5s
Running query:   RUNNING 264.2s
Running query:   RUNNING 264.9s
Running query:   RUNNING 265.6s
Running query:   RUNNING 266.3s
Running query:   RUNNING 267.4s
Running query:   RUNNING 268.5s
Running query:   RUNNING 271.0s
Running query:   RUNNING 272.0s
Running query:   RUNNING 273.1s
Running query:   RUNNING 273.9s
Running query:   RUNNING 274.8s
Running query:   RUNNING 275.5s
Running query:   RUNNING 276.3s
Running query:   RUNNING 277.1s
Running query:   RUNNING 278.3s
Running query:   RUNNING 279.1s
Running query:   RUNNING 279.8s
Running query:   RUNNING 280.5s
Running query:   RUNNING 281.2s
Running query:   RUNNING 283.4s
Running query:   RUNNING 284.1s
Running query:   RUNNING 284.8s
Running query:   RUNNING 285.5s
Running query:   RUNNING 286.3s
Running query:   RUNNING 287.0s
Running query:   RUNNING 288.4s
Running query:   RUNNING 289.4s
Running query:   RUNNING 290.1s
Running query:   RUNNING 290.9s
Running query:   RUNNING 291.6s
Running query:   RUNNING 292.3s
Running query:   RUNNING 293.1s
Running query:   RUNNING 293.8s
Running query:   RUNNING 294.5s
Running query:   RUNNING 295.3s
Running query:   RUNNING 296.0s
Running query:   RUNNING 297.2s
Running query:   RUNNING 298.2s
Running query:   RUNNING 298.9s
Running query:   RUNNING 299.6s
Running query:   RUNNING 300.4s
Running query:   RUNNING 301.1s
Running query:   RUNNING 302.0s
Running query:   RUNNING 302.7s
Running query:   RUNNING 303.4s
Running query:   RUNNING 304.6s
Running query:   RUNNING 305.3s
Running query:   RUNNING 306.1s
Running query:   RUNNING 307.5s
Running query:   RUNNING 308.4s
Running query:   RUNNING 309.7s
Running query:   RUNNING 310.4s
Running query:   RUNNING 311.1s
Running query:   RUNNING 312.0s
Running query:   RUNNING 313.0s
Running query:   RUNNING 315.4s
Running query:   RUNNING 316.6s
Running query:   RUNNING 317.9s
Running query:   RUNNING 318.6s
Running query:   RUNNING 319.3s
Running query:   RUNNING 320.0s
Running query:   RUNNING 320.7s
Running query:   RUNNING 321.4s
Running query:   RUNNING 322.1s
Running query:   RUNNING 322.9s
Running query:   RUNNING 323.6s
Running query:   RUNNING 324.3s
Running query:   RUNNING 325.0s
Running query:   RUNNING 325.7s
Running query:   RUNNING 326.4s
Running query:   RUNNING 327.2s
Running query:   RUNNING 327.9s
Running query:   RUNNING 328.6s
Running query:   RUNNING 329.4s
Running query:   RUNNING 330.2s
Running query:   RUNNING 330.9s
Running query:   RUNNING 331.6s
Running query:   RUNNING 332.3s
Running query:   RUNNING 333.0s
Running query:   RUNNING 333.8s
Running query:   RUNNING 334.5s
Running query:   RUNNING 335.2s
Running query:   RUNNING 335.9s
Running query:   RUNNING 336.6s
Running query:   RUNNING 337.5s
Running query:   RUNNING 338.2s
Running query:   RUNNING 338.9s
Running query:   RUNNING 339.6s
Running query:   RUNNING 340.3s
Running query:   RUNNING 341.0s
Running query:   RUNNING 341.7s
Running query:   RUNNING 342.5s
Running query:   RUNNING 343.2s
Running query:   RUNNING 344.0s
Running query:   RUNNING 345.4s
Running query:   RUNNING 346.1s
Running query:   RUNNING 347.2s
Running query:   RUNNING 348.0s
Running query:   RUNNING 349.2s
Running query:   RUNNING 350.7s
Running query:   RUNNING 351.4s
Running query:   RUNNING 352.5s
Running query:   RUNNING 353.2s
Running query:   RUNNING 354.5s
Running query:   RUNNING 355.6s
Running query:   RUNNING 356.3s
Running query:   RUNNING 357.0s
Running query:   RUNNING 357.8s
Running query:   RUNNING 358.5s
Running query:   RUNNING 359.2s
Running query:   RUNNING 359.9s
Running query:   RUNNING 360.6s
Running query:   RUNNING 361.3s
Running query:   RUNNING 362.1s
Running query:   RUNNING 362.8s
Running query:   RUNNING 363.5s
Running query:   RUNNING 364.2s
Running query:   RUNNING 365.0s
Running query:   RUNNING 365.7s
Running query:   RUNNING 366.7s
Running query:   RUNNING 367.4s
Running query:   RUNNING 368.1s
Running query:   RUNNING 368.9s
Running query:   RUNNING 369.6s
Running query:   RUNNING 370.4s
Running query:   RUNNING 371.1s
Running query:   RUNNING 371.8s
Running query:   RUNNING 372.6s
Running query:   RUNNING 373.4s
Running query:   RUNNING 374.1s
Running query:   RUNNING 374.9s
Running query:   RUNNING 375.6s
Running query:   RUNNING 376.3s
Running query:   RUNNING 377.0s
Running query:   RUNNING 377.8s
Running query:   RUNNING 378.5s
Running query:   RUNNING 379.2s
Running query:   RUNNING 380.0s
Running query:   RUNNING 380.7s
Running query:   RUNNING 381.4s
Running query:   RUNNING 382.3s
Running query:   RUNNING 383.3s
Running query:   RUNNING 384.0s
Running query:   RUNNING 384.8s
Running query:   RUNNING 385.9s
Running query:   RUNNING 387.2s
Running query:   RUNNING 388.8s
Running query:   RUNNING 390.4s
Running query:   RUNNING 391.1s
Running query:   RUNNING 391.9s
Running query:   RUNNING 393.1s
Running query:   RUNNING 393.8s
Running query:   RUNNING 394.5s
Running query:   RUNNING 395.2s
Running query:   RUNNING 395.9s
Running query:   RUNNING 396.7s
Running query:   RUNNING 397.4s
Running query:   RUNNING 398.1s
Running query:   RUNNING 398.9s
Running query:   RUNNING 399.9s
Running query:   RUNNING 400.7s
Running query:   RUNNING 401.4s
Running query:   RUNNING 402.2s
Running query:   RUNNING 403.0s
Running query:   RUNNING 403.7s
Running query:   RUNNING 404.4s
Running query:   RUNNING 405.2s
Running query:   RUNNING 405.9s
Running query:   RUNNING 406.6s
```

```r
plate = substr(concordanceResult$sample_id, 1, 9)
concordanceResult = cbind(concordanceResult, plate)
```


```r
ggplot(concordanceResult) +
  geom_point(aes(x=sample_id, y=concordance, color=plate), size=4) +
  xlab("Sample") +
  ylab("Concordance") +
  ggtitle("Concordance with Genotyping Data") +
  scale_colour_brewer(name="Sample Prep Plate", palette="Set1") +
  theme(axis.text.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        title=element_text(size=24),
        legend.text=element_text(size=20)) 
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
  scale_colour_brewer(name="Sample Prep Plate", palette="Set1") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=22),
        title=element_text(size=24),
        legend.text=element_text(size=20)) 
```

<img src="figure/pca-batch-effect-publication-1.png" title="plot of chunk pca-batch-effect-publication" alt="plot of chunk pca-batch-effect-publication" style="display: block; margin: auto;" />

#### IBS

```r
ibsDataFlowFilename = '/Users/gmcinnes/data/all-genomes-ibs-2.tsv'
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
    theme_bw() +
    theme(axis.ticks=element_blank(), 
          axis.text=element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          title=element_text(size=24)) +
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
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender), size=4) +
  xlab("Sample") +
  ylab("Heterozygosity Rate") +
  ggtitle("Heterozygosity Rate on the X Chromosome") +
  scale_colour_brewer(palette="Set1", name="Gender") +
  theme(axis.text.x=if(nrow(result) <= 20)
    {element_text(angle = 90, hjust = 1)} else {element_blank()},
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=20),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title=element_text(size=22),
    title=element_text(size=24),
    legend.text=element_text(size=20)) 
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
  stat_smooth() +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  scale_x_continuous(expand = c(0, 0), labels=comma) +
  ggtitle("Ti/Tv by 100,000 base pair windows\non Chromosome 1") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        title=element_text(size=24))
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

Retrieving data: 19.5s
Retrieving data: 27.0s
Retrieving data: 40.3s
Retrieving data: 66.3s
Retrieving data: 89.5s
Retrieving data: 100.0s
Retrieving data: 139.5s
Retrieving data: 170.2s
Retrieving data: 177.6s
```


```r
ggplot(titv, aes(x=average_depth, y=titv_ratio, color=call_call_set_name)) + 
  geom_point(size=2) +
  ggtitle("Ti/Tv By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        title=element_text(size=24))
```

<img src="figure/titv-by-depth-publication-1.png" title="plot of chunk titv-by-depth-publication" alt="plot of chunk titv-by-depth-publication" style="display: block; margin: auto;" />


