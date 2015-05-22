SELECT
  CHROM,
  POS,
  ID,
  REF,
  ALT,
  ROUND(QUAL, 2) AS QUAL,
  "PASS" AS FILTER,
  "." AS INFO,
  sample_calls,
FROM js(
  (SELECT
     reference_name, 
     start, 
     reference_bases, 
     alternate_bases, 
     names, 
     call.call_set_name, 
     call.genotype, 
     call.DP, 
     call.FILTER, 
     call.QUAL
   FROM
     [va_aaa_pilot_data.all_genomes_expanded_vcfs_java3]
   OMIT RECORD IF EVERY(call.genotype <= 0) # OR COUNT(names) = 0
   #LIMIT 1000
  ),
  reference_name, 
  start, 
  reference_bases, 
  alternate_bases, 
  names, 
  call.call_set_name, 
  call.genotype, 
  call.DP, 
  call.FILTER, 
  call.QUAL,
  "[{name: 'debug', type: 'string'},
    {name: 'CHROM', type: 'string'},
    {name: 'POS', type: 'integer'},
    {name: 'ID', type: 'string'},
    {name: 'REF', type: 'string'},
    {name: 'ALT', type: 'string'},
    {name: 'QUAL', type: 'float'},

    {name: 'sample_calls', type: 'string'}]",
#    {name: 'FILTER', type: 'string'},    
#    {name: 'INFO', type: 'string'},
  "function(r, emit) {
     var alt = r.alternate_bases.join();
     if(1 == alt.length && 1 == r.reference_bases.length) { 
       var samples = '';
       var quality = 0.0;
       var filter = '';
       for (call in r.call) {
         samples += r.call[call].call_set_name + ':' + r.call[call].genotype.join('/') + '|' ;
         quality += r.call[call].QUAL;
       }
       emit({
         debug: JSON.stringify(r),
         CHROM: r.reference_name,
         POS: r.start,
         ID: r.names.join([separator = '|']),
         REF: r.reference_bases,
         ALT: alt,
         QUAL: quality/r.call.length,
         sample_calls: samples
       });
     }
   }")