# Publication Plots - QC

## Setup


Plot theme

```r
plot_theme = theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(axis.line = element_line(colour = "black"))
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
Running query:   RUNNING  2.0s
Running query:   RUNNING  2.7s
Running query:   RUNNING  3.4s
Running query:   RUNNING  4.0s
Running query:   RUNNING  4.7s
Running query:   RUNNING  5.4s
Running query:   RUNNING  6.3s
Running query:   RUNNING  7.1s
Running query:   RUNNING  7.9s
Running query:   RUNNING  8.6s
Running query:   RUNNING  9.3s
Running query:   RUNNING 10.3s
Running query:   RUNNING 11.6s
Running query:   RUNNING 12.8s
Running query:   RUNNING 13.5s
Running query:   RUNNING 14.4s
Running query:   RUNNING 15.3s
Running query:   RUNNING 16.0s
Running query:   RUNNING 17.2s
Running query:   RUNNING 18.0s
Running query:   RUNNING 18.7s
Running query:   RUNNING 19.5s
Running query:   RUNNING 20.8s
Running query:   RUNNING 21.5s
Running query:   RUNNING 22.2s
Running query:   RUNNING 22.9s
Running query:   RUNNING 23.6s
Running query:   RUNNING 24.4s
Running query:   RUNNING 25.1s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.6s
Running query:   RUNNING 27.4s
Running query:   RUNNING 28.1s
Running query:   RUNNING 28.9s
Running query:   RUNNING 29.6s
Running query:   RUNNING 30.4s
Running query:   RUNNING 31.2s
Running query:   RUNNING 31.9s
Running query:   RUNNING 32.7s
Running query:   RUNNING 33.4s
Running query:   RUNNING 34.2s
Running query:   RUNNING 34.9s
Running query:   RUNNING 35.7s
Running query:   RUNNING 36.4s
Running query:   RUNNING 37.2s
Running query:   RUNNING 37.9s
Running query:   RUNNING 38.7s
Running query:   RUNNING 39.5s
Running query:   RUNNING 40.2s
Running query:   RUNNING 40.9s
Running query:   RUNNING 41.7s
Running query:   RUNNING 42.5s
Running query:   RUNNING 44.0s
Running query:   RUNNING 45.1s
Running query:   RUNNING 46.0s
Running query:   RUNNING 47.2s
Running query:   RUNNING 48.0s
Running query:   RUNNING 48.7s
Running query:   RUNNING 49.6s
Running query:   RUNNING 50.3s
Running query:   RUNNING 51.1s
Running query:   RUNNING 51.9s
Running query:   RUNNING 52.6s
Running query:   RUNNING 53.4s
Running query:   RUNNING 54.7s
Running query:   RUNNING 55.4s
Running query:   RUNNING 56.2s
Running query:   RUNNING 56.9s
Running query:   RUNNING 57.7s
Running query:   RUNNING 58.3s
Running query:   RUNNING 59.0s
Running query:   RUNNING 59.7s
Running query:   RUNNING 60.5s
Running query:   RUNNING 61.2s
Running query:   RUNNING 62.0s
Running query:   RUNNING 62.7s
Running query:   RUNNING 63.5s
Running query:   RUNNING 64.2s
Running query:   RUNNING 65.0s
Running query:   RUNNING 65.7s
Running query:   RUNNING 66.4s
Running query:   RUNNING 67.2s
Running query:   RUNNING 67.9s
Running query:   RUNNING 68.6s
Running query:   RUNNING 69.4s
Running query:   RUNNING 73.6s
Running query:   RUNNING 74.4s
Running query:   RUNNING 75.1s
Running query:   RUNNING 75.8s
Running query:   RUNNING 76.6s
Running query:   RUNNING 77.8s
Running query:   RUNNING 78.8s
Running query:   RUNNING 79.6s
Running query:   RUNNING 92.5s
Running query:   RUNNING 99.6s
Running query:   RUNNING 101.5s
Running query:   RUNNING 102.3s
Running query:   RUNNING 103.0s
Running query:   RUNNING 104.3s
Running query:   RUNNING 105.1s
Running query:   RUNNING 106.5s
Running query:   RUNNING 107.4s
Running query:   RUNNING 108.5s
Running query:   RUNNING 109.9s
Running query:   RUNNING 111.1s
Running query:   RUNNING 112.6s
Running query:   RUNNING 113.4s
Running query:   RUNNING 114.0s
Running query:   RUNNING 114.7s
Running query:   RUNNING 115.5s
Running query:   RUNNING 116.2s
Running query:   RUNNING 117.0s
Running query:   RUNNING 117.7s
Running query:   RUNNING 118.5s
Running query:   RUNNING 119.2s
Running query:   RUNNING 119.9s
Running query:   RUNNING 120.7s
Running query:   RUNNING 121.4s
Running query:   RUNNING 122.2s
Running query:   RUNNING 122.9s
Running query:   RUNNING 123.6s
Running query:   RUNNING 124.4s
Running query:   RUNNING 125.1s
Running query:   RUNNING 125.9s
Running query:   RUNNING 126.6s
Running query:   RUNNING 127.3s
Running query:   RUNNING 128.1s
Running query:   RUNNING 128.8s
Running query:   RUNNING 129.6s
Running query:   RUNNING 130.5s
Running query:   RUNNING 131.2s
Running query:   RUNNING 132.0s
Running query:   RUNNING 132.8s
Running query:   RUNNING 133.5s
Running query:   RUNNING 134.2s
Running query:   RUNNING 134.9s
Running query:   RUNNING 135.7s
Running query:   RUNNING 136.4s
Running query:   RUNNING 137.2s
Running query:   RUNNING 138.0s
Running query:   RUNNING 139.0s
Running query:   RUNNING 139.9s
Running query:   RUNNING 140.7s
Running query:   RUNNING 141.8s
Running query:   RUNNING 142.5s
Running query:   RUNNING 143.2s
Running query:   RUNNING 144.1s
Running query:   RUNNING 144.8s
Running query:   RUNNING 146.1s
Running query:   RUNNING 146.9s
Running query:   RUNNING 148.7s
Running query:   RUNNING 149.5s
Running query:   RUNNING 150.2s
Running query:   RUNNING 150.8s
Running query:   RUNNING 151.6s
Running query:   RUNNING 152.3s
Running query:   RUNNING 153.1s
Running query:   RUNNING 153.8s
Running query:   RUNNING 154.6s
Running query:   RUNNING 155.3s
Running query:   RUNNING 156.0s
Running query:   RUNNING 156.8s
Running query:   RUNNING 157.5s
Running query:   RUNNING 158.2s
Running query:   RUNNING 159.0s
Running query:   RUNNING 159.7s
Running query:   RUNNING 160.4s
Running query:   RUNNING 161.2s
Running query:   RUNNING 161.9s
Running query:   RUNNING 162.6s
Running query:   RUNNING 163.4s
Running query:   RUNNING 164.1s
Running query:   RUNNING 164.8s
Running query:   RUNNING 165.6s
Running query:   RUNNING 166.3s
Running query:   RUNNING 167.0s
Running query:   RUNNING 167.7s
Running query:   RUNNING 168.5s
Running query:   RUNNING 169.2s
Running query:   RUNNING 169.9s
Running query:   RUNNING 170.7s
Running query:   RUNNING 171.4s
Running query:   RUNNING 172.1s
Running query:   RUNNING 172.9s
Running query:   RUNNING 173.6s
Running query:   RUNNING 174.4s
Running query:   RUNNING 175.1s
Running query:   RUNNING 175.9s
Running query:   RUNNING 176.6s
Running query:   RUNNING 177.4s
Running query:   RUNNING 178.1s
Running query:   RUNNING 178.9s
Running query:   RUNNING 179.6s
Running query:   RUNNING 180.4s
Running query:   RUNNING 181.3s
Running query:   RUNNING 182.1s
Running query:   RUNNING 182.9s
Running query:   RUNNING 183.7s
Running query:   RUNNING 184.5s
Running query:   RUNNING 185.3s
Running query:   RUNNING 186.0s
Running query:   RUNNING 187.3s
Running query:   RUNNING 188.0s
Running query:   RUNNING 188.7s
Running query:   RUNNING 189.4s
Running query:   RUNNING 190.1s
Running query:   RUNNING 190.9s
Running query:   RUNNING 191.7s
Running query:   RUNNING 192.4s
Running query:   RUNNING 193.1s
Running query:   RUNNING 193.9s
Running query:   RUNNING 194.6s
Running query:   RUNNING 195.3s
Running query:   RUNNING 196.1s
Running query:   RUNNING 196.8s
Running query:   RUNNING 197.5s
Running query:   RUNNING 198.3s
Running query:   RUNNING 199.0s
Running query:   RUNNING 199.8s
Running query:   RUNNING 200.5s
Running query:   RUNNING 201.2s
Running query:   RUNNING 202.0s
Running query:   RUNNING 202.8s
Running query:   RUNNING 203.5s
Running query:   RUNNING 204.2s
Running query:   RUNNING 205.0s
Running query:   RUNNING 205.7s
Running query:   RUNNING 207.2s
Running query:   RUNNING 208.0s
Running query:   RUNNING 208.7s
Running query:   RUNNING 209.4s
Running query:   RUNNING 210.1s
Running query:   RUNNING 210.9s
Running query:   RUNNING 211.6s
Running query:   RUNNING 212.3s
Running query:   RUNNING 213.1s
Running query:   RUNNING 214.4s
Running query:   RUNNING 215.1s
Running query:   RUNNING 215.9s
Running query:   RUNNING 216.9s
Running query:   RUNNING 217.6s
Running query:   RUNNING 218.4s
Running query:   RUNNING 219.2s
Running query:   RUNNING 219.9s
Running query:   RUNNING 220.7s
Running query:   RUNNING 221.4s
Running query:   RUNNING 223.9s
Running query:   RUNNING 225.1s
Running query:   RUNNING 225.7s
Running query:   RUNNING 226.4s
Running query:   RUNNING 227.1s
Running query:   RUNNING 227.8s
Running query:   RUNNING 228.5s
Running query:   RUNNING 229.3s
Running query:   RUNNING 230.0s
Running query:   RUNNING 230.7s
Running query:   RUNNING 231.5s
Running query:   RUNNING 232.2s
Running query:   RUNNING 232.9s
Running query:   RUNNING 233.6s
Running query:   RUNNING 234.4s
Running query:   RUNNING 235.1s
Running query:   RUNNING 235.8s
Running query:   RUNNING 236.6s
Running query:   RUNNING 237.3s
Running query:   RUNNING 238.0s
Running query:   RUNNING 238.8s
Running query:   RUNNING 239.5s
Running query:   RUNNING 240.2s
Running query:   RUNNING 240.9s
Running query:   RUNNING 241.7s
Running query:   RUNNING 242.4s
Running query:   RUNNING 243.1s
Running query:   RUNNING 243.8s
Running query:   RUNNING 244.6s
Running query:   RUNNING 245.3s
Running query:   RUNNING 246.1s
Running query:   RUNNING 246.8s
Running query:   RUNNING 247.5s
Running query:   RUNNING 248.3s
Running query:   RUNNING 249.0s
Running query:   RUNNING 249.8s
Running query:   RUNNING 250.5s
Running query:   RUNNING 251.3s
Running query:   RUNNING 252.4s
Running query:   RUNNING 253.6s
Running query:   RUNNING 255.0s
Running query:   RUNNING 258.0s
Running query:   RUNNING 258.8s
Running query:   RUNNING 259.6s
Running query:   RUNNING 263.3s
Running query:   RUNNING 263.9s
Running query:   RUNNING 264.6s
Running query:   RUNNING 265.4s
Running query:   RUNNING 266.1s
Running query:   RUNNING 266.8s
Running query:   RUNNING 267.5s
Running query:   RUNNING 268.3s
Running query:   RUNNING 269.0s
Running query:   RUNNING 269.7s
Running query:   RUNNING 270.5s
Running query:   RUNNING 271.2s
Running query:   RUNNING 271.9s
Running query:   RUNNING 272.7s
Running query:   RUNNING 273.4s
Running query:   RUNNING 274.1s
Running query:   RUNNING 274.9s
Running query:   RUNNING 275.6s
Running query:   RUNNING 276.3s
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
  scale_colour_brewer(name="Sample Prep Plate", palette="Set1") +
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
  scale_colour_brewer(name="Sample Prep Plate", palette="Set1") +
  plot_theme +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position = c(0.85,0.8))
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

Retrieving data:  6.1s
Retrieving data: 11.8s
Retrieving data: 18.9s
Retrieving data: 32.3s
Retrieving data: 39.8s
Retrieving data: 42.6s
Retrieving data: 44.9s
Retrieving data: 48.4s
Retrieving data: 52.6s
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


