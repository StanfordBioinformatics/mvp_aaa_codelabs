
class Queries(object):
    GENDER_CHECK = "gender-check-fail.sql"
    GENOTYPING_CONCORDANCE = "genotyping-conconcordance-fail.sql"
    INBREEDING_COEFFICIENT = "inbreeding-coefficient-fail.sql"
    INBREEDING_COEFFICIENT_METRICS = "inbreeding-coefficient-metrics.sql"
    MISSINGNESS_SAMPLE_LEVEL = "missingness-sample-level-fail.sql"
    PRIVATE_VARIANTS = "private-variants-fail.sql"
    PRIVATE_VARIANT_METRICS = "private-variants-metrics.sql"
    VARIANT_DEPTH = "variant-depth-fail.sql"
    MISSINGNESS_VARIANT_LEVEL = "variant-level-missingness-fail.sql"

    SAMPLE_LEVEL_QC_QUERIES = [GENDER_CHECK,
                               #GENOTYPING_CONCORDANCE,
                               INBREEDING_COEFFICIENT,
                               MISSINGNESS_SAMPLE_LEVEL,
                               PRIVATE_VARIANTS
    ]

    VARIANT_LEVEL_QC_QUERIES = [VARIANT_DEPTH,
                                MISSINGNESS_VARIANT_LEVEL]

    AVERAGE_STDDEV = {
        PRIVATE_VARIANTS: PRIVATE_VARIANT_METRICS,
        INBREEDING_COEFFICIENT: INBREEDING_COEFFICIENT_METRICS
    }

    PRESET_CUTOFFS = {
        GENOTYPING_CONCORDANCE: {"_CUTOFF_": "0.99"},
        GENDER_CHECK: {"_MALE_CUTOFF_": "0.2",
                       "_FEMALE_CUTOFF_": "0.5"},
        MISSINGNESS_SAMPLE_LEVEL: {"_CUTOFF_": "1"},
        MISSINGNESS_VARIANT_LEVEL: {"_CUTOFF_": "0.9"},
        VARIANT_DEPTH: {"_MAX_VALUE_": "70",
                        "_MIN_VALUE_": "8"}
    }