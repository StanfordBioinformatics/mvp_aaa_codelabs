# mvp_aaa_codelabs
This codelab goes over the procedures used to prepare and analyze a pilot set of genomes from the Million Veterans Program on Google Cloud.  The goal of this project is to develop a robust, scalable solution for performing large-scale genomic analysis.  We have taken analyses that are traditionally performed on an HPC and implemented them on Google Cloud so that they run on distributed file systems.  The result is a set of tools that run exceptionally fast when compared to standard genomic analysis tools and can be scaled to a large number of genomes.

The work presented here is featured in a [manuscript](https://academic.oup.com/bioinformatics/article/4036385/Cloud-based-Interactive-Analytics-for-Terabytes-of?searchresult=1) we published in Bioinformatics.

### Part 1: [BigQuery Setup](./BigQuery-Setup.md)
This part of the codelab goes over the necessary steps to create tables in BigQuery with genomic information.  

### Part 2: [Sample Level Quality Control](./Sample-Level-QC.md)
This section explains and demonstrates the methods used to perform quality control on each sample in our cohort.

### Part 3: [Variant Level Quality Control](./Variant-Level-QC.md) 
This section explains and demonstrates the methods used to perfrom quality control on each variant in the genomes.

###  Part 4: [QC Implementation](./QC-Implementation.md)
This part of the codelab demonstrates execution of the QC pipeline, removing low quality samples, and flagging low quality variants.


The overall workflow discussed here is visualized in the workflow below.

<img src="figure/bioinformatics_workflow.png" title="bioinformatics workflow" alt="bioinformatics workflow" style="display: block; margin: auto;" />

