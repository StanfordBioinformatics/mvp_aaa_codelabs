# Get all variants that have a chi squared value above a definited limit
SELECT
  variant_id,
  chisq,
  "hardy_weinberg" AS failure_reason,
FROM js(
    (SELECT
      variant_id,
      reference_name,
      start,
      reference_bases,
      alternate_bases,
      hom_ref AS obs_hom1,
      het AS obs_het,
      hom_alt AS obs_hom2,
      hom_ref + het + hom_alt AS sample_count,
    FROM (
      SELECT
        variant_id,
        reference_name,
        start,
        END,
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        SUM(EVERY(0 = call.genotype)) WITHIN call AS hom_ref,
        SUM(EVERY(1 = call.genotype)) WITHIN call AS hom_alt,
        SUM(SOME(0 = call.genotype)
          AND SOME(1 = call.genotype)) WITHIN call AS het,
      FROM
        [_THE_EXPANDED_TABLE_]
      # Optionally add a clause here to limit the query to a particular
      # region of the genome.
      #_WHERE_
      HAVING
        # Skip 1/2 genotypes
        num_alts = 1
        )),
        // Start javascript function
        // Input Columns
        variant_id, reference_name, start, reference_bases, alternate_bases, obs_hom1, obs_hom2, obs_het, sample_count,
        // Output Schema
        "[{name: 'variant_id', type: 'string'},
        {name: 'reference_name', type: 'string'},
        {name: 'start', type: 'integer'},
        {name: 'reference_bases', type: 'string'},
        {name: 'alternate_bases', type: 'string'},
        {name: 'obs_hom1', type: 'integer'},
        {name: 'obs_het', type: 'integer'},
        {name: 'obs_hom2', type: 'integer'},
        {name: 'e_hom1', type: 'integer'},
        {name: 'e_het', type: 'integer'},
        {name: 'e_hom2', type: 'integer'},        
        {name: 'chisq', type: 'float'}]",
        // Function
        "function(r, emit) {
          var e_hom1 = Math.pow((r.obs_hom1 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
          var e_het = 2 * ((r.obs_hom1 + (r.obs_het/2)) / r.sample_count) * ((r.obs_hom2 + (r.obs_het/2)) / r.sample_count) * r.sample_count;
          var e_hom2 = Math.pow((r.obs_hom2 + (r.obs_het/2)) / r.sample_count, 2) * r.sample_count;
          var chisq = (Math.pow(r.obs_hom1 - e_hom1, 2) / e_hom1) + (Math.pow(r.obs_het - e_het, 2) / e_het) + (Math.pow(r.obs_hom2 - e_hom2, 2) / e_hom2);
          emit({
            variant_id: r.variant_id,
            reference_name: r.reference_name,
            start: r.start,
            reference_bases: r.reference_bases,
            alternate_bases: r.alternate_bases,
            obs_hom1: r.obs_hom1,
            obs_hom2: r.obs_hom2,
            obs_het: r.obs_het,
            e_hom1: e_hom1,
            e_hom2: e_hom2,
            e_het: e_het,
            chisq: chisq
          })
        }"
      )
WHERE
  chisq > _CUTOFF_
#_ORDER_BY_

