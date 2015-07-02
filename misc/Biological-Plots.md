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
                          "_GENOTYPING_TABLE_"="va_aaa_pilot_data.genotyping_data",
                          "_ANNOVAR_TABLE_"="resources.annovar_hg19_aaa",
                          "_CHR_LENGTHS_"="stanford.edu:gbsc-stanford-google:resources.hg19_Assembly_BinaRuns")
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
counts = ggplot(genotypeCountResult) +
  geom_bar(aes(x=Genotype, y=Cnt), stat="identity") +
  xlab("Genotypes") + 
  ylab("SNV Count") + 
  scale_y_continuous(labels=comma, expand = c(0, 0)) + 
  plot_theme
counts
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
snvs = ggplot(data=snps, aes(y=Cnt, x=reference_name)) + 
  geom_point(size=4) + 
  ylab("SNV Count") +
  xlab("Chromosome") +
  scale_y_continuous(label=comma) +
  plot_theme +
  theme(axis.text.x = element_text(angle=90, vjust=1)) 
snvs
```

<img src="figure/snv-counts-publication-1.png" title="plot of chunk snv-counts-publication" alt="plot of chunk snv-counts-publication" style="display: block; margin: auto;" />


```r
indels = ggplot(data=indels, aes(y=Cnt, x=reference_name)) + 
  geom_point(size=4) + 
  ylab("Indel Count") +
  xlab("Chromosome") +
  scale_y_continuous(label=comma) +
  plot_theme +
  theme(axis.text.x = element_text(angle=90, vjust=1)) 
indels
```

<img src="figure/indel-counts-publication-1.png" title="plot of chunk indel-counts-publication" alt="plot of chunk indel-counts-publication" style="display: block; margin: auto;" />

#### Saturation Rate

```r
# Need to redo queries for this
genomeCount = c(1,2,3,4,5,10,50,100,200,300,400,478)
snvCount = c(3590360,4847512,5627244,6158953,6616457,
              8014799,11841547,14387937,18693833,21567571,23638061,25890797)
saturationRate = data_frame(genomeCount, snvCount)
```



```r
saturation = ggplot(saturationRate) +
  geom_point(aes(x=genomeCount, y=snvCount), size=4) +
  xlab("Number of Genomes") +
  ylab("Unique SNVs") +
  scale_y_continuous(label=comma) +
  plot_theme
saturation
```

<img src="figure/saturation-publication-1.png" title="plot of chunk saturation-publication" alt="plot of chunk saturation-publication" style="display: block; margin: auto;" />

## Multiplot

```r
multiplot(snvs,indels, counts, saturation, cols=2)
```

<img src="figure/variant-multiplot-publication-1.png" title="plot of chunk variant-multiplot-publication" alt="plot of chunk variant-multiplot-publication" style="display: block; margin: auto;" />

#### Ti/Tv for SNVs

```r
titvBySample <- DisplayAndDispatchQuery("../sql/titv-by-sample.sql",
                                              project=project,
                                              replacements=queryReplacements)
```

```
SELECT
call.call_set_name AS sample.id,
transitions,
transversions,
transitions/transversions AS titv,
num_variants,
FROM 
(
  SELECT
  call.call_set_name,
  SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
  SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                   'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
  COUNT(mutation) AS num_variants
  FROM 
  (
    SELECT
    call.call_set_name,
    CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
    COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    GROUP_CONCAT(QC) WITHIN RECORD AS qc,
    GROUP_CONCAT(call.QC) WITHIN CALL AS call_qc,
    FROM
    FLATTEN([va_aaa_pilot_data.multi_sample_variants_full_qc], alternate_bases)
    OMIT call IF EVERY (call.genotype <= 0)
    OR SOME(call.qc IS NOT NULL)
    HAVING
    num_alts = 1
    AND QC IS NULL
  )
  GROUP BY
  call.call_set_name
)
```


```r
titv = ggplot(titvBySample) +
  geom_boxplot(aes(x="SNV", y=titv)) +
  xlab("Sample") +
  ylab("ti/tv") +
  theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(axis.text.x=element_blank())
titv
```

<img src="figure/titv-by-sample-1.png" title="plot of chunk titv-by-sample" alt="plot of chunk titv-by-sample" style="display: block; margin: auto;" />

#### Het/Hom for SNVs

```r
hetHomSnv <- DisplayAndDispatchQuery("../sql/hethom-snv.sql",
                                              project=project,
                                              replacements=queryReplacements)
```

```
SELECT
sample.id,
(hom_AA_count + het_RA_count) AS all_SNV_count,
hom_AA_count,
het_RA_count,
ROUND((het_RA_count/hom_AA_count)*1000)/1000 AS Het_Hom_ratio
FROM
(
  SELECT
  call.call_set_name AS sample.id,
  SUM(genotype="1/1" AND (LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")))) AS hom_AA_count,
  SUM((genotype="0/1" OR genotype="1/0" OR genotype="1/2") AND (LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")))) AS het_RA_count    
  FROM
  (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(STRING(call.genotype), "/") WITHIN call AS genotype,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    call.FILTER
    FROM 
    FLATTEN([va_aaa_pilot_data.multi_sample_variants_full_qc], alternate_bases)
    OMIT call IF SOME(call.qc IS NOT NULL)
    HAVING
    cohort_qc IS NULL AND call.FILTER = "PASS"
  ) 
  GROUP BY
  sample.id
)
GROUP BY
sample.id,
hom_AA_count,
het_RA_count,
all_SNV_count,
Het_Hom_ratio
ORDER BY
sample.id;
```


```r
hethom= ggplot(hetHomSnv) +
  geom_boxplot(aes(x="SNV", y=Het_Hom_ratio)) +
  xlab("Sample") +
  ylab("het/hom") +
  theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(axis.text.x=element_blank())
hethom
```

<img src="figure/hethom-snv-1.png" title="plot of chunk hethom-snv" alt="plot of chunk hethom-snv" style="display: block; margin: auto;" />

#### SNVs in dbSNP

```r
dbSNPcount <- DisplayAndDispatchQuery("../sql/dbsnp-count.sql",
                                              project=project,
                                              replacements=queryReplacements)
```

```
SELECT
call.call_set_name AS sample_id,
VAR_type,
COUNT(call.call_set_name) AS num_VAR,
COUNT(dbSNP) AS num_VAR_dbSNP
FROM (
  SELECT
  call.call_set_name,
  IF(LENGTH(reference_bases)=1 AND (LENGTH(alternates)=1 OR (LENGTH(alternates)=3 AND alternates CONTAINS ",")) , "SNV", "INDEL") AS VAR_type,
  dbSNP
  FROM (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(names) WITHIN record AS dbSNP,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    call.FILTER
    FROM 
    FLATTEN([va_aaa_pilot_data.multi_sample_variants_full_qc], alternate_bases)
    OMIT call IF SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL AND call.FILTER = "PASS"
  )
)
GROUP BY
sample_id,
VAR_type
ORDER BY
sample_id,
VAR_type
```

```r
dbSNPcount$proportion = dbSNPcount$num_VAR_dbSNP/dbSNPcount$num_VAR
```


```r
dbSNP = ggplot(dbSNPcount) +
  geom_boxplot(aes(x=VAR_type, y=proportion)) +
  xlab("Sample") +
  ylab("Proportion of Variants in dbSNP") +
  scale_y_continuous(label=comma) +
  theme_minimal(base_size = 18, base_family = "Helvetica") 
dbSNP
```

<img src="figure/dbsnp-snv-1.png" title="plot of chunk dbsnp-snv" alt="plot of chunk dbsnp-snv" style="display: block; margin: auto;" />

#### Private SNV Counts

```r
privateSNVs <- DisplayAndDispatchQuery("../sql/private-snv-counts.sql",
                                              project=project,
                                              replacements=queryReplacements)
```

```
SELECT
call.call_set_name AS sample_id,
COUNT(call.call_set_name) AS private_SNVs_count
FROM
(
  SELECT
  call.call_set_name,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples_with_VARs,
  GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
  GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
  FROM
  [va_aaa_pilot_data.multi_sample_variants_full_qc]
  OMIT 
  call IF EVERY(call.FILTER != "PASS") 
  OR SOME(call.QC IS NOT NULL)
  HAVING
  num_samples_with_VARs = 1
  AND cohort_qc IS NULL
)
GROUP EACH BY
sample_id
ORDER BY
sample_id ASC;
```


```r
privateVariants = ggplot(privateSNVs) +
  geom_boxplot(aes(x="SNV", y=private_SNVs_count)) +
  xlab("Sample") +
  ylab("Private Variants") +
  theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(axis.text.x=element_blank())
privateVariants 
```

<img src="figure/private-snv-1.png" title="plot of chunk private-snv" alt="plot of chunk private-snv" style="display: block; margin: auto;" />

## Multiplot

```r
multiplot(titv, dbSNP, hethom, privateVariants, cols=2)
```

<img src="figure/boxplots-multiplot-publication-1.png" title="plot of chunk boxplots-multiplot-publication" alt="plot of chunk boxplots-multiplot-publication" style="display: block; margin: auto;" />

#### Variant Rarity

```r
rarity <- DisplayAndDispatchQuery("../sql/variant-rarity.sql",
                                              project=project,
                                              replacements=queryReplacements)
```

```
SELECT 
reference_name,
rarity,
COUNT(rarity) AS count
FROM (
  SELECT 
  reference_name, 
  start,
  end,
  reference_bases,
  alternates,
  num_samples,
  ROUND(num_samples/459, 4) AS allele_frequency,
  CASE WHEN num_samples = 1 THEN "very_rare"
  WHEN num_samples = 2 THEN "rare" 
  WHEN num_samples >= 3 AND num_samples <= 23 THEN "relatively_common"
  WHEN num_samples >= 24 THEN "common"
  ELSE "NONE"
  END AS rarity 
  FROM (
    SELECT
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc,
    COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
    FROM
    [va_aaa_pilot_data.multi_sample_variants_full_qc]
    OMIT 
    call IF EVERY(call.FILTER != "PASS")
    OR SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL 
    AND reference_name in ("chr1", "chr2", "chr3", "chr4", "chr5",
                           "chr6", "chr7", "chr8", "chr9", "chr10",
                           "chr11", "chr12", "chr13", "chr14", "chr15",
                           "chr16", "chr17", "chr18", "chr19", "chr20",
                           "chr21", "chr22", "chrX", "chrY")
  )
  GROUP EACH BY
  reference_name,
  start,
  end,
  reference_bases,
  alternates,
  call_qc,
  cohort_qc,
  num_samples,
  allele_frequency,
  rarity
)
GROUP BY 
reference_name, 
rarity
ORDER BY
reference_name,
rarity
```

```r
rarity$total = sapply(1:nrow(rarity), function(x){
  chr = rarity[x,]$reference_name
  variantCount = sum(rarity[rarity$reference_name == chr,]$count)
  variantCount
})
rarity$proportion = rarity$count/rarity$total
rarity$reference_name <- factor(rarity$reference_name, levels=chromosomes)
rarity$rarity <- factor(rarity$rarity, levels = c("very_rare", "rare", "relatively_common", "common"), ordered=TRUE)
```


```r
ggplot(rarity) +
  geom_bar(aes(factor(reference_name), proportion, fill=rarity, order=rarity), stat='identity', position='stack') +
  guides(fill = guide_legend(reverse=TRUE)) +
  xlab("Chromosome") +
  ylab("Proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette=6) +
  theme_minimal(base_size = 18, base_family = "Helvetica") + 
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        legend.title=element_blank())
```

<img src="figure/variant-rarity-1.png" title="plot of chunk variant-rarity" alt="plot of chunk variant-rarity" style="display: block; margin: auto;" />

```r
rarityPlot
```

```
Error in data.frame(x = c("chr12", "chr9", "chrX", "chr9", "chr6", "chr19", : arguments imply differing number of rows: 96, 0
```

#### Regional Variant Count

```r
regionalCounts <- DisplayAndDispatchQuery("../sql/regional-variant-count.sql",
                                              project=project,
                                              replacements=queryReplacements,
                                              max=Inf)
```

```
SELECT
Sample_id,
Chromosome,
region,
cnt,
LENGTH AS chromosome_length
FROM
(
  SELECT
  call.call_set_name AS Sample_id,
  reference_name AS Chromosome,
  Func_refGene AS region,
  COUNT(geno.start) AS cnt
  FROM 
  [resources.annovar_hg19_aaa] AS anno
  JOIN EACH
  (
    SELECT
    call.call_set_name,
    reference_name,
    start,
    end,
    reference_bases,
    GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternates,
    GROUP_CONCAT(call.QC) WITHIN call AS call_qc,
    GROUP_CONCAT(QC) WITHIN RECORD AS cohort_qc
    FROM
    [va_aaa_pilot_data.multi_sample_variants_full_qc]
    OMIT 
    call IF EVERY(call.FILTER != "PASS") 
    OR SOME(call.QC IS NOT NULL)
    HAVING
    cohort_qc IS NULL 
  ) AS geno
  ON
  anno.Chr = geno.reference_name
  AND anno.Start = geno.start
  AND anno.End = geno.end
  AND anno.ALT = geno.alternates
  GROUP EACH BY
  Sample_id,
  Chromosome,
  region
) AS cal
JOIN
[stanford.edu:gbsc-stanford-google:resources.hg19_Assembly_BinaRuns] AS ref
ON
cal.Chromosome = ref.CHR
ORDER BY
Sample_id,
region,
Chromosome;
Retrieving data:  3.2sRetrieving data:  4.6sRetrieving data:  6.3sRetrieving data:  7.6sRetrieving data:  8.9sRetrieving data: 10.3sRetrieving data: 11.6sRetrieving data: 13.2sRetrieving data: 14.6sRetrieving data: 15.8sRetrieving data: 17.0sRetrieving data: 18.2sRetrieving data: 19.6sRetrieving data: 21.1sRetrieving data: 22.8sRetrieving data: 24.1s
```

```r
regionalCounts$Chromosome <- factor(regionalCounts$Chromosome, levels=chromosomes)
regionalCounts$scaled <- regionalCounts$cnt/regionalCounts$chromosome_length
regionalCounts <- regionalCounts[complete.cases(regionalCounts),]
```


```r
exonic = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="exonic"), aes(Chromosome,scaled)) +
  ggtitle("Variants in Exonic Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
exonic
```

<img src="figure/exonic-1.png" title="plot of chunk exonic" alt="plot of chunk exonic" style="display: block; margin: auto;" />


```r
intronic = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="intronic"), aes(Chromosome,scaled)) +
  ggtitle("Variants in Intronic Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
intronic
```

<img src="figure/intronic-1.png" title="plot of chunk intronic" alt="plot of chunk intronic" style="display: block; margin: auto;" />


```r
utr3 = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="UTR3"), aes(Chromosome,scaled)) +
  ggtitle("Variants in UTR3 Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
utr3
```

<img src="figure/utr3-1.png" title="plot of chunk utr3" alt="plot of chunk utr3" style="display: block; margin: auto;" />


```r
utr5 = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="UTR5"), aes(Chromosome,scaled)) +
  ggtitle("Variants in UTR5 Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
utr5
```

<img src="figure/utr5-1.png" title="plot of chunk utr5" alt="plot of chunk utr5" style="display: block; margin: auto;" />


```r
intergenic = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="intergenic"), aes(Chromosome,scaled)) +
  ggtitle("Variants in Intergenic Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
intergenic
```

<img src="figure/intergenic-1.png" title="plot of chunk intergenic" alt="plot of chunk intergenic" style="display: block; margin: auto;" />


```r
splicing = ggplot(regionalCounts) +
  geom_boxplot(data=subset(regionalCounts,region=="splicing"), aes(Chromosome,scaled)) +
  ggtitle("Variants in Splicing Regions") +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank())
splicing
```

<img src="figure/splicing-1.png" title="plot of chunk splicing" alt="plot of chunk splicing" style="display: block; margin: auto;" />

#### Multiplot

```r
multiplot(exonic,utr3, intergenic, intronic, utr5, splicing, cols=2)
```

<img src="figure/regional-counts-1.png" title="plot of chunk regional-counts" alt="plot of chunk regional-counts" style="display: block; margin: auto;" />
