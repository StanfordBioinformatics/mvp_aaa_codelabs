# Missingness Data Preparation on BigQuery

This document describers how data was prepared in BigQuery to check sample missingness.

## Create a table with all reference calls

A BigQuery [query](./sql/select_ref_calls_brca1.sql) was run to generate a table that lists every position within BRCA1 for each of the 5 sample genomes as either a reference call or not.

This table was exported to a table called 5_genomes_ref_calls_brca1 within the qc_tables dataset.

## Create a table with all variant calls.

A BigQuery [query](./sql/select_variants_brca1.sql) was run to generate a table that lists every variant within BRCA1 for each of the 5 sample genomes.

This table was exported to a table called 5_genomes_variants_brca1 within the qc_tables dataset.