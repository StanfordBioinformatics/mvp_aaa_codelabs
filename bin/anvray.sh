#!/bin/bash
date
hostname
module load annovar/20141112
vcf=$1
root=`echo $vcf | rev | cut -f1 -d/ | rev | cut -f1 -d.`

table_annovar.pl $vcf $ANNOVAR/humandb/ -buildver hg19 -out $root.annovar -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gerp++elem,gerp++gt2,esp6500si_all,ALL.sites.2014_10,snp138,ljb26_all,cg69,cosmic70,clinvar_20140929,nci60,exac02,caddgt10,avsift,popfreq_all -operation g,g,g,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
