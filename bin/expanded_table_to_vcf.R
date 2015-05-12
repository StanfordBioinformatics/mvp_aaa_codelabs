# Generate pseudo multi sample vcfs for each chromosome from an expanded variant table
setwd("/Users/gmcinnes/GitHub/mvp_aaa_codelabs/qc")
project <- "gbsc-gcp-project-mvp"
source("./rHelpers/setup.R")
tableReplacement <- list("_THE_EXPANDED_TABLE_"="va_aaa_pilot_data.all_genomes_expanded_vcfs_java2")
for (c in c(seq(1,22),"X","Y")) {
  chr = paste("'chr", c, "'", sep="")
  chrSub = list("_CHR_"=chr)
  outputTable = paste("qc_tables.vcf_chr", c, sep="")
  print(chr)
  print(outputTable)
  result <- DisplayAndDispatchQuery("./sql/multisample-vcf.sql",
                                  project=project,
                                  replacements=c(tableReplacement,chrSub),
                                  outputTable=outputTable)
}
