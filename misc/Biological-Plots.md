# Publication Plots - Biological

## Setup


Plot theme

```r
plot_theme = theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(axis.line = element_line(colour = "black"))
```

## Biological Query Plots
Tables for biological queries

```r
queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.genome_calls_full_qc",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.multi_sample_variants_full_qc",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions",
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data")
```



#### Genotype Counts

```r
genotypeCountResult <- DisplayAndDispatchQuery("../sql/genotype-counts.sql",
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
    GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotype,
    GROUP_CONCAT(QC) WITHIN RECORD AS qc,
    GROUP_CONCAT(call.QC) WITHIN CALL AS call_qc,
    FROM 
    [va_aaa_pilot_data.multi_sample_variants_full_qc]
    OMIT call IF SOME(call.qc IS NOT NULL)
    HAVING QC IS NULL
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
  scale_y_continuous(labels=comma, expand = c(0, 0)) + 
  plot_theme
```

<img src="figure/genotype-counts-publication-1.png" title="plot of chunk genotype-counts-publication" alt="plot of chunk genotype-counts-publication" style="display: block; margin: auto;" />

#### Variant Counts By Chromosome

```r
variantCountResult <- DisplayAndDispatchQuery("../sql/variants-by-chromosome.sql",
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
    GROUP_CONCAT(QC) WITHIN RECORD AS qc,
    GROUP_CONCAT(call.QC) WITHIN CALL AS call_qc,
    IF(LENGTH(reference_bases)=1 AND LENGTH(alternate_bases)=1, "SNV", "INDEL") AS VAR_type
    FROM 
    [va_aaa_pilot_data.multi_sample_variants_full_qc]
    OMIT call IF EVERY(call.genotype <= 0) 
      OR SOME(call.qc IS NOT NULL)
    HAVING QC IS NULL
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

Running query:   RUNNING  2.4s
Running query:   RUNNING  3.0s
Running query:   RUNNING  3.6s
Running query:   RUNNING  4.2s
Running query:   RUNNING  4.8s
Running query:   RUNNING  5.4s
Running query:   RUNNING  6.1s
Running query:   RUNNING  6.7s
Running query:   RUNNING  7.3s
Running query:   RUNNING  7.9s
Running query:   RUNNING  8.5s
Running query:   RUNNING  9.1s
Running query:   RUNNING  9.7s
Running query:   RUNNING 10.4s
Running query:   RUNNING 11.0s
Running query:   RUNNING 11.6s
Running query:   RUNNING 12.2s
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
ggplot(data=snps, aes(y=Cnt, x=reference_name)) + 
  geom_point(size=4) + 
  ylab("SNV Count") +
  xlab("Chromosome") +
  scale_y_continuous(label=comma) +
  plot_theme +
  theme(axis.text.x = element_text(angle=90, vjust=1)) 
```

<img src="figure/snv-counts-publication-1.png" title="plot of chunk snv-counts-publication" alt="plot of chunk snv-counts-publication" style="display: block; margin: auto;" />


```r
ggplot(data=indels, aes(y=Cnt, x=reference_name)) + 
  geom_point(size=4) + 
  ylab("Indel Count") +
  xlab("Chromosome") +
  scale_y_continuous(label=comma) +
  plot_theme +
  theme(axis.text.x = element_text(angle=90, vjust=1)) 
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
  plot_theme
```

<img src="figure/saturation-publication-1.png" title="plot of chunk saturation-publication" alt="plot of chunk saturation-publication" style="display: block; margin: auto;" />
