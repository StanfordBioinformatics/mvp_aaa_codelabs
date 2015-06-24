# Poster Plots

setwd("/Users/gmcinnes/GitHub/mvp_aaa_codelabs/qc")

# Set the Google Cloud Platform project id under which these queries will run.
project <- "gbsc-gcp-project-mvp"

# Set up for BigQuery access.
source("./rHelpers/setup.R")

queryReplacements <- list("_THE_TABLE_"="va_aaa_pilot_data.all_genomes_gvcfs",
                          "_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java2",
                          "_BLACKLISTED_TABLE_"="resources.blacklisted_positions")
sampleData <- read.csv("./data/patient_info.csv")
sampleInfo <- select(sampleData, call_call_set_name=Catalog.ID, gender=Gender)

# Ti/Tv By Depth
query <- "./sql/ti-tv-by-depth.sql"
titv <- DisplayAndDispatchQuery(query,
                                  project=project,
                                  replacements=c(tableReplacement))

ggplot(titv, aes(x=average_depth, y=titv_ratio, color=call_call_set_name), alpha=0.01) + 
  geom_point() +
  ggtitle("Ti/Tv Ratio By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        title=element_text(size=24,face="bold"))

#my_breaks = c(5, 50, 200, 350, 500)
#my_breaks = c(100, 200, 300, 400, 500)

#ggplot(titv, aes(x=average_depth, y=titv_ratio)) + 
  #geom_hex() +
  #stat_binhex(bins=100) +
  #ggtitle("Ti/Tv Ratio By Depth") +
  #xlab("Coverage Depth") + 
  #ylab("Ti/Tv") +
  #theme_bw() +
  #scale_fill_gradientn(colours=c("red","blue"), name = "count", trans = "log",
  #                     breaks = my_breaks, labels = my_breaks) +
  #theme(axis.line = element_line(colour = "black"),
  #      panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
  #      panel.border = element_blank(),
  #      panel.background = element_blank(),
  #      legend.position = "none",
  #      axis.text=element_text(size=16),
  #      axis.title=element_text(size=20),
  #      title=element_text(size=24,face="bold"))


### IBS 
require(reshape2)
require(dplyr)
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

RemoveSamples <- function(ibsData) {
  toRemove = c('LP6005692-DNA_F12',
               'LP6005243-DNA_G12',
               'LP6005243-DNA_A12',
               'LP6005243-DNA_G11',
               'LP6005243-DNA_C12',
               'LP6005243-DNA_D12',
               'LP6005243-DNA_F12',
               'LP6005243-DNA_H12',
               'LP6005693-DNA_B03',
               'LP6005243-DNA_E12',
               'LP6005243-DNA_H11',
               'LP6005692-DNA_G12',
               'LP6005144-DNA_E04',
               'LP6005243-DNA_F11',
               'LP6005243-DNA_B12',
               'LP6005693-DNA_C03',
               'LP6005692-DNA_H12',
               'LP6005693-DNA_D03',
               'LP6005693-DNA_H01',
               'LP6005695-DNA_A01')
  ibsData <- ibsData[!ibsData$sample1 %in% toRemove,]
  ibsData <- ibsData[!ibsData$sample2 %in% toRemove,]
  return(ibsData)
}

ibsDataflowDataSample <- RemoveSamples(ibsDataflowDataSample)

SampleIBSMatrix <- function(ibsData, sampleSize=50) {
  individuals <- unique(ibsData$sample1)
  sample <- sample(individuals, sampleSize)
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
DrawHeatMap(ibsDataflowDataSubset)



titvWindowResults <- DisplayAndDispatchQuery("./sql/ti-tv-ratio.sql",
                                  project=project,
                                  replacements=c("#_WHERE_"="WHERE reference_name = 'chr1'",
                                                 "_WINDOW_SIZE_"="100000",
                                                 queryReplacements))


ggplot(titvWindowResults, aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous() +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by 100,000 base pair windows on Chromosome 1")


#### Sex inference
result <- DisplayAndDispatchQuery("./sql/gender-check.sql",
                                  project=project,
                                  replacements=tableReplacement)

joinedResult <- inner_join(result, sampleInfo)


ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender), size=5) +
  xlab("Sample") +
  ylab("Heterozygosity Rate ") +
  ggtitle("Heterozygosity Rate on the X Chromosome") +
  scale_colour_brewer(palette="Set1") +
  theme(axis.text.x=if(nrow(result) <= 20)
    {element_text(angle = 90, hjust = 1)} else {element_blank()},
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=16),
    axis.title=element_text(size=20),
    title=element_text(size=24,face="bold")) 

####

plate = substr(concordanceResult$sample_id, 1, 9)
concordanceResult = cbind(concordanceResult, plate)

# COLORED
ggplot(concordanceResult) +
  geom_point(aes(x=sample_id, y=concordance, color=plate), size=4) +
  xlab("Sample") +
  ylab("Concordance") +
  ggtitle("Concordance with Genotyping Data") +
  scale_colour_brewer(name="Sample Preparation Plate", palette="Set1") +
  theme(axis.text.x=element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=24),
    title=element_text(size=28,face="bold")) 

# BLACK
ggplot(concordanceResult) +
  geom_point(aes(x=sample_id, y=concordance), size=4) +
  xlab("Sample") +
  ylab("Concordance") +
  ggtitle("Concordance with Genotyping Data") +
  scale_colour_brewer(name="Sample Preparation Plate", palette="Set1") +
  theme(axis.text.x=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        title=element_text(size=28,face="bold")) 



pcaFile = '/Users/gmcinnes/data/pca-all-genomes-all-references-no-tissue.tsv'
pcaResult = read.table(pcaFile)
names(pcaResult) = c('sample_id','pc1', 'pc2', 'something')
plate = substr(pcaResult$sample_id, 1, 9)
pcaResult = cbind(pcaResult, plate)

ggplot(pcaResult, aes(pc1, pc2, color=plate)) + 
  geom_point(size=7) +
  ggtitle("Principal Component Analysis") +
  xlab("Principal Component 1") + 
  ylab("Principal Component 2") +
  scale_colour_brewer(name="Sample Preparation Plate", palette="Set1") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        title=element_text(size=28,face="bold"),
        legend.text=element_text(size=20)) 

#########################################

## CONCORDANCE
# COLORED
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
        title=element_text(size=32,face="bold"),
        legend.text=element_text(size=24)) 

# PCA
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

# TiTv
ggplot(titv, aes(x=average_depth, y=titv_ratio, color=call_call_set_name)) + 
  geom_point(size=3) +
  ggtitle("Ti/Tv Ratio By Depth") +
  xlab("Coverage Depth") + 
  ylab("Ti/Tv") +
  #theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        title=element_text(size=32,face="bold"))

# Sex
ggplot(joinedResult) +
  geom_point(aes(x=call_call_set_name, y=perct_het_alt_in_snvs, color=gender), size=5) +
  xlab("Sample") +
  ylab("Heterozygosity Rate ") +
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
  axis.title=element_text(size=28),
  title=element_text(size=32,face="bold"),
  legend.text=element_text(size=24)) 

# Ti/Tv By Genomic Window
require(scales)
ggplot(titvWindowResults, aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
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


# IBS
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
DrawHeatMap(ibsDataflowDataSubset)
