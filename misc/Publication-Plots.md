# Publication Plots

## Setup

```r
require(scales)
require(reshape2)
require(dplyr)
require(plyr)

# Set Working Directory
setwd("/Users/gmcinnes/src/mvp_aaa_codelabs/qc")

# Set the Google Cloud Platform project id under which these queries will run.
project <- "gbsc-gcp-project-mvp"

# Set up for BigQuery access.
source("./rHelpers/setup.R")

# Define Query Substitutions
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java2",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions",
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data")

# Read in sample information
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

chromosomeLengths <- read.table("http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes")
names(chromosomeLengths) = c("reference_name", "length")
```


## Biological Query Plots
#### Genotype Counts

```r
genotypeCountResult <- DisplayAndDispatchQuery("./sql/genotype-counts.sql",
                                             project=project,
                                             replacements=queryReplacements)
```

```
SELECT
Genotype,
COUNT(genotype) AS Cnt
FROM
(
  SELECT
  reference_name,
  start,
  reference_bases,
  alternates,
  genotype
  FROM
  (
    SELECT
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotype
    FROM 
    [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
  )
  GROUP EACH BY
  reference_name,
  start,
  reference_bases,
  alternates,
  genotype
)
GROUP BY
Genotype
ORDER BY
Genotype
```

```r
genotypeCountResult <- genotypeCountResult[complete.cases(genotypeCountResult),]
```


```r
ggplot(genotypeCountResult) +
  geom_bar(aes(x=Genotype, y=Cnt), stat="identity") +
  xlab("Genotypes") + 
  ylab("SNV Count") + 
  ggtitle("SNV Counts for Each Genotype") +
  scale_y_continuous(labels=comma, expand = c(0, 0)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
```

<img src="figure/genotype-counts-publication-1.png" title="plot of chunk genotype-counts-publication" alt="plot of chunk genotype-counts-publication" style="display: block; margin: auto;" />

#### Variant Counts By Chromosome

```r
variantCountResult <- DisplayAndDispatchQuery("./sql/variants-by-chromosome.sql",
                                             project=project,
                                             replacements=queryReplacements)
```

```
SELECT
reference_name,
VAR_type,
COUNT(VAR_type) AS Cnt
FROM
(
  SELECT
  reference_name,
  start,
  reference_bases,
  alternates,
  VAR_type
  FROM
  (
    SELECT
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    IF(LENGTH(reference_bases)=1 AND LENGTH(alternate_bases)=1, "SNV", "INDEL") AS VAR_type
    FROM 
    [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
    OMIT call IF EVERY(call.genotype <= 0)
  )
  GROUP EACH BY
  reference_name,
  start,
  reference_bases,
  alternates,
  VAR_type
)
GROUP BY
reference_name,
VAR_type
ORDER BY
reference_name,
VAR_type
```

```r
variantCountResult <- join(variantCountResult, chromosomeLengths, by = "reference_name")
variantCountResult$scaled_count <- variantCountResult$Cnt / variantCountResult$length

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
snps <- variantCountResult[grep("SNV", variantCountResult$VAR_type), ]
snps$reference_name <- factor(snps$reference_name, levels=chromosomes)
snps <- snps[complete.cases(snps),]

indels <- variantCountResult[grep("INDEL", variantCountResult$VAR_type), ]
indels$reference_name <- factor(indels$reference_name, levels=chromosomes)
indels <- indels[complete.cases(indels),]
```


```r
ggplot(data=snps, aes(y=scaled_count, x=reference_name)) + 
  geom_point() + 
  ylab("Scaled SNV Count") +
  xlab("Chromosome") +
  scale_y_continuous(label=comma) +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(angle=90, vjust=1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
```

<img src="figure/snv-counts-publication-1.png" title="plot of chunk snv-counts-publication" alt="plot of chunk snv-counts-publication" style="display: block; margin: auto;" />


```r
ggplot(data=indels, aes(y=scaled_count, x=reference_name)) + 
  geom_point() + 
  ylab("Scaled Indel Count") +
  xlab("chromosome") +
  scale_y_continuous(label=comma) +
  theme(text = element_text(size=12), 
        axis.text.x = element_text(angle=90, vjust=1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
```

<img src="figure/indel-counts-publication-1.png" title="plot of chunk indel-counts-publication" alt="plot of chunk indel-counts-publication" style="display: block; margin: auto;" />

#### Saturation Rate

```r
# Need to redo queries for this
genome_count = c(1,2,3,4,5,10,50,100,200,300,400,478)
snv_count = c(3590360,4847512,5627244,6158953,6616457,
              8014799,11841547,14387937,18693833,21567571,23638061,25890797)
saturation_rate = data_frame(genome_count, snv_count)
```

```r
ggplot(saturation_rate) +
  geom_point(aes(x=genome_count, y=snv_count), size=4) +
  xlab("Number of Genomes") +
  ylab("Unique SNVs") +
  scale_y_continuous(label=comma) +
  ggtitle("Saturation Rate") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
```

<img src="figure/saturation-publication-1.png" title="plot of chunk saturation-publication" alt="plot of chunk saturation-publication" style="display: block; margin: auto;" />

## QC Plots
#### Genotyping Concordance

```r
concordanceResult <- DisplayAndDispatchQuery("./sql/genotyping-concordance.sql",
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
          [va_aaa_pilot_data.all_genomes_gvcfs]
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
Running query:   RUNNING  6.3s
Running query:   RUNNING  6.9s
Running query:   RUNNING  7.5s
Running query:   RUNNING  8.1s
Running query:   RUNNING  8.9s
Running query:   RUNNING  9.5s
Running query:   RUNNING 10.2s
Running query:   RUNNING 10.8s
Running query:   RUNNING 11.4s
Running query:   RUNNING 12.0s
Running query:   RUNNING 12.7s
Running query:   RUNNING 13.3s
Running query:   RUNNING 13.9s
Running query:   RUNNING 14.5s
Running query:   RUNNING 15.1s
Running query:   RUNNING 15.8s
Running query:   RUNNING 16.4s
Running query:   RUNNING 17.2s
Running query:   RUNNING 17.8s
Running query:   RUNNING 18.4s
Running query:   RUNNING 19.0s
Running query:   RUNNING 19.7s
Running query:   RUNNING 20.3s
Running query:   RUNNING 20.9s
Running query:   RUNNING 21.5s
Running query:   RUNNING 22.2s
Running query:   RUNNING 22.8s
Running query:   RUNNING 23.4s
Running query:   RUNNING 24.0s
Running query:   RUNNING 24.6s
Running query:   RUNNING 25.3s
Running query:   RUNNING 25.9s
Running query:   RUNNING 26.5s
Running query:   RUNNING 27.2s
Running query:   RUNNING 27.8s
Running query:   RUNNING 28.4s
Running query:   RUNNING 29.2s
Running query:   RUNNING 29.8s
Running query:   RUNNING 30.5s
Running query:   RUNNING 31.1s
Running query:   RUNNING 31.8s
Running query:   RUNNING 32.4s
Running query:   RUNNING 33.0s
Running query:   RUNNING 33.7s
Running query:   RUNNING 34.3s
Running query:   RUNNING 34.9s
Running query:   RUNNING 35.6s
Running query:   RUNNING 36.2s
Running query:   RUNNING 36.8s
Running query:   RUNNING 37.4s
Running query:   RUNNING 38.0s
Running query:   RUNNING 38.7s
Running query:   RUNNING 39.3s
Running query:   RUNNING 39.9s
Running query:   RUNNING 40.5s
Running query:   RUNNING 41.2s
Running query:   RUNNING 41.8s
Running query:   RUNNING 42.5s
Running query:   RUNNING 43.1s
Running query:   RUNNING 43.7s
Running query:   RUNNING 44.3s
Running query:   RUNNING 44.9s
Running query:   RUNNING 45.6s
Running query:   RUNNING 46.2s
Running query:   RUNNING 46.8s
Running query:   RUNNING 47.5s
Running query:   RUNNING 48.1s
Running query:   RUNNING 48.7s
Running query:   RUNNING 49.3s
Running query:   RUNNING 50.0s
Running query:   RUNNING 50.6s
Running query:   RUNNING 51.2s
Running query:   RUNNING 51.8s
Running query:   RUNNING 52.5s
Running query:   RUNNING 53.1s
Running query:   RUNNING 53.7s
Running query:   RUNNING 54.3s
Running query:   RUNNING 54.9s
Running query:   RUNNING 55.6s
Running query:   RUNNING 56.2s
Running query:   RUNNING 56.8s
Running query:   RUNNING 57.5s
Running query:   RUNNING 58.1s
Running query:   RUNNING 58.7s
Running query:   RUNNING 59.4s
Running query:   RUNNING 60.0s
Running query:   RUNNING 60.6s
Running query:   RUNNING 61.2s
Running query:   RUNNING 61.8s
Running query:   RUNNING 62.5s
Running query:   RUNNING 63.1s
Running query:   RUNNING 63.8s
Running query:   RUNNING 64.4s
Running query:   RUNNING 65.0s
Running query:   RUNNING 65.6s
Running query:   RUNNING 66.3s
Running query:   RUNNING 66.9s
Running query:   RUNNING 67.5s
Running query:   RUNNING 68.1s
Running query:   RUNNING 68.8s
Running query:   RUNNING 69.4s
Running query:   RUNNING 70.0s
Running query:   RUNNING 70.6s
Running query:   RUNNING 71.2s
Running query:   RUNNING 71.9s
Running query:   RUNNING 72.5s
Running query:   RUNNING 73.1s
Running query:   RUNNING 73.8s
Running query:   RUNNING 74.4s
Running query:   RUNNING 75.0s
Running query:   RUNNING 75.7s
Running query:   RUNNING 76.3s
Running query:   RUNNING 76.9s
Running query:   RUNNING 77.5s
Running query:   RUNNING 78.2s
Running query:   RUNNING 78.8s
Running query:   RUNNING 79.4s
Running query:   RUNNING 80.0s
Running query:   RUNNING 80.6s
Running query:   RUNNING 81.3s
Running query:   RUNNING 81.9s
Running query:   RUNNING 82.5s
Running query:   RUNNING 83.1s
Running query:   RUNNING 83.8s
Running query:   RUNNING 84.4s
Running query:   RUNNING 85.0s
Running query:   RUNNING 85.6s
Running query:   RUNNING 86.3s
Running query:   RUNNING 86.9s
Running query:   RUNNING 87.6s
Running query:   RUNNING 88.2s
Running query:   RUNNING 88.8s
Running query:   RUNNING 89.4s
Running query:   RUNNING 90.0s
Running query:   RUNNING 90.7s
Running query:   RUNNING 91.3s
Running query:   RUNNING 91.9s
Running query:   RUNNING 92.5s
Running query:   RUNNING 93.1s
Running query:   RUNNING 93.8s
Running query:   RUNNING 94.4s
Running query:   RUNNING 95.0s
Running query:   RUNNING 95.6s
Running query:   RUNNING 96.2s
Running query:   RUNNING 96.9s
Running query:   RUNNING 97.5s
Running query:   RUNNING 98.2s
Running query:   RUNNING 98.8s
Running query:   RUNNING 99.5s
Running query:   RUNNING 100.1s
Running query:   RUNNING 100.8s
Running query:   RUNNING 101.4s
Running query:   RUNNING 102.1s
Running query:   RUNNING 102.7s
Running query:   RUNNING 103.4s
Running query:   RUNNING 104.0s
Running query:   RUNNING 104.7s
Running query:   RUNNING 105.3s
Running query:   RUNNING 106.0s
Running query:   RUNNING 106.6s
Running query:   RUNNING 107.3s
Running query:   RUNNING 107.9s
Running query:   RUNNING 108.6s
Running query:   RUNNING 109.2s
Running query:   RUNNING 109.8s
Running query:   RUNNING 110.5s
Running query:   RUNNING 111.1s
Running query:   RUNNING 111.7s
Running query:   RUNNING 112.3s
Running query:   RUNNING 112.9s
Running query:   RUNNING 113.5s
Running query:   RUNNING 114.2s
Running query:   RUNNING 114.8s
Running query:   RUNNING 115.4s
Running query:   RUNNING 116.0s
Running query:   RUNNING 116.7s
Running query:   RUNNING 117.4s
Running query:   RUNNING 118.0s
Running query:   RUNNING 118.6s
Running query:   RUNNING 119.4s
Running query:   RUNNING 120.0s
Running query:   RUNNING 120.7s
Running query:   RUNNING 121.3s
Running query:   RUNNING 122.2s
Running query:   RUNNING 122.8s
Running query:   RUNNING 123.4s
Running query:   RUNNING 124.0s
Running query:   RUNNING 124.7s
Running query:   RUNNING 125.3s
Running query:   RUNNING 125.9s
Running query:   RUNNING 126.7s
Running query:   RUNNING 127.3s
Running query:   RUNNING 128.0s
Running query:   RUNNING 128.6s
Running query:   RUNNING 129.2s
Running query:   RUNNING 129.8s
Running query:   RUNNING 130.4s
Running query:   RUNNING 131.1s
Running query:   RUNNING 131.7s
Running query:   RUNNING 132.3s
Running query:   RUNNING 132.9s
Running query:   RUNNING 133.6s
Running query:   RUNNING 134.2s
Running query:   RUNNING 134.8s
Running query:   RUNNING 135.4s
Running query:   RUNNING 136.1s
Running query:   RUNNING 136.7s
Running query:   RUNNING 137.3s
Running query:   RUNNING 137.9s
Running query:   RUNNING 138.6s
Running query:   RUNNING 139.2s
Running query:   RUNNING 139.8s
Running query:   RUNNING 140.5s
Running query:   RUNNING 141.1s
Running query:   RUNNING 141.8s
Running query:   RUNNING 142.4s
Running query:   RUNNING 143.1s
Running query:   RUNNING 143.7s
Running query:   RUNNING 144.3s
Running query:   RUNNING 145.0s
Running query:   RUNNING 145.6s
Running query:   RUNNING 146.3s
Running query:   RUNNING 146.9s
Running query:   RUNNING 147.6s
Running query:   RUNNING 148.2s
Running query:   RUNNING 148.8s
Running query:   RUNNING 149.5s
Running query:   RUNNING 150.2s
Running query:   RUNNING 150.8s
Running query:   RUNNING 151.4s
Running query:   RUNNING 152.1s
Running query:   RUNNING 152.7s
Running query:   RUNNING 153.4s
Running query:   RUNNING 154.1s
Running query:   RUNNING 154.7s
Running query:   RUNNING 155.3s
Running query:   RUNNING 156.0s
Running query:   RUNNING 156.6s
Running query:   RUNNING 157.3s
Running query:   RUNNING 157.9s
Running query:   RUNNING 158.6s
Running query:   RUNNING 159.2s
Running query:   RUNNING 159.8s
Running query:   RUNNING 160.5s
Running query:   RUNNING 161.1s
Running query:   RUNNING 161.8s
Running query:   RUNNING 162.4s
Running query:   RUNNING 163.1s
Running query:   RUNNING 163.7s
Running query:   RUNNING 164.3s
Running query:   RUNNING 165.0s
Running query:   RUNNING 165.6s
Running query:   RUNNING 166.3s
Running query:   RUNNING 166.9s
Running query:   RUNNING 167.6s
Running query:   RUNNING 168.2s
Running query:   RUNNING 168.8s
Running query:   RUNNING 169.5s
Running query:   RUNNING 170.1s
Running query:   RUNNING 170.8s
Running query:   RUNNING 171.4s
Running query:   RUNNING 172.1s
Running query:   RUNNING 172.7s
Running query:   RUNNING 173.4s
Running query:   RUNNING 174.0s
Running query:   RUNNING 174.7s
Running query:   RUNNING 175.3s
Running query:   RUNNING 176.0s
Running query:   RUNNING 176.6s
Running query:   RUNNING 177.3s
Running query:   RUNNING 177.9s
Running query:   RUNNING 178.6s
Running query:   RUNNING 179.2s
Running query:   RUNNING 179.9s
Running query:   RUNNING 180.5s
Running query:   RUNNING 181.2s
Running query:   RUNNING 181.8s
Running query:   RUNNING 182.5s
Running query:   RUNNING 183.2s
Running query:   RUNNING 183.8s
Running query:   RUNNING 184.5s
Running query:   RUNNING 185.1s
Running query:   RUNNING 185.8s
Running query:   RUNNING 186.4s
Running query:   RUNNING 187.1s
Running query:   RUNNING 187.7s
Running query:   RUNNING 188.4s
Running query:   RUNNING 189.0s
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
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
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
  geom_point(size=7) +
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
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 
```

<img src="figure/pca-batch-effect-publication-1.png" title="plot of chunk pca-batch-effect-publication" alt="plot of chunk pca-batch-effect-publication" style="display: block; margin: auto;" />

#### Ti/Tv By Depth

```r
query <- "./sql/ti-tv-by-depth.sql"
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
        [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
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

Retrieving data:  2.1s
Retrieving data:  3.3s
Retrieving data:  4.6s
Retrieving data:  5.8s
Retrieving data:  7.0s
Retrieving data:  7.9s
Retrieving data:  8.9s
Retrieving data:  9.8s
```


```r
ggplot(titv, aes(x=average_depth, y=titv_ratio, color=call_call_set_name)) + 
  geom_point(size=3) +
  ggtitle("Ti/Tv Ratio By Depth") +
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
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"))
```

<img src="figure/titv-by-depth-publication-1.png" title="plot of chunk titv-by-depth-publication" alt="plot of chunk titv-by-depth-publication" style="display: block; margin: auto;" />

#### Sex Inference

```r
result <- DisplayAndDispatchQuery("./sql/gender-check.sql",
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
      [va_aaa_pilot_data.all_genomes_expanded_vcfs_java2]
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

Running query:   RUNNING  2.5s
Running query:   RUNNING  3.1s
```

```r
joinedResult <- inner_join(result, sampleInfo)
```


```r
ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender), size=5) +
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
    axis.text=element_text(size=24),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title=element_text(size=28),
    title=element_text(size=32,face="bold"),
    legend.text=element_text(size=24)) 
```

<img src="figure/sex-inference-publication-1.png" title="plot of chunk sex-inference-publication" alt="plot of chunk sex-inference-publication" style="display: block; margin: auto;" />

#### Ti/Tv By Genomic Window

```r
titvWindowResults <- DisplayAndDispatchQuery("./sql/ti-tv-ratio.sql",
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
      [va_aaa_pilot_data.all_genomes_gvcfs]
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
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle("Ti/Tv by 100,000 base pair windows\non Chromosome 1") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"))
```

<img src="figure/titv-by-genomic-window-publication-1.png" title="plot of chunk titv-by-genomic-window-publication" alt="plot of chunk titv-by-genomic-window-publication" style="display: block; margin: auto;" />

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
          title=element_text(size=24,face="bold")) +
    geom_tile(aes(fill=ibsScore), colour="white") +
    scale_fill_gradient(low="white", high="steelblue", na.value="black",
                        guide=guide_colourbar(title= "IBS Score")) +
    labs(list(title="Identity By State (IBS) Heat Map",
              x="Sample", y="Sample")) 
  p 
}
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
          axis.title=element_text(size=28),
          legend.text=element_text(size=24),
          title=element_text(size=32,face="bold")) +
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
