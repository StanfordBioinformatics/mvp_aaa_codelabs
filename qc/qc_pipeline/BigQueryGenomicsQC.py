import os
import sys
import re
from GenomicsQueries import Queries
from BigQueryClient import BigQuery
from config import Config
from GoogleGenomicsClient import GoogleGenomicsClient
import logging

class GenomicsQC(object):
    def __init__(self, verbose=False):
        self.query_repo = Config.QUERY_REPO
        self.variant_table = Config.VARIANT_TABLE
        self.expanded_table = Config.EXPANDED_TABLE
        self.client_secrets_path = Config.CLIENT_SECRETS
        self.project_number = Config.PROJECT_NUMBER
        self.setup_log(verbose)
        self.bq = BigQuery()

    #### Specific types of QC functions ####
    # Sample level QC
    def sample_qc(self):
        logging.info("Running Sample Level QC")
        queries = Queries.SAMPLE_LEVEL_QC_QUERIES
        for q in queries:
            logging.debug(q)
            self.run_analysis(query_file=q, level='sample')

    # Variant level QC
    def variant_qc(self):
        logging.info("Running Variant Level QC")
        queries = Queries.VARIANT_LEVEL_QC_QUERIES
        for q in queries:
            logging.debug(q)
            self.run_analysis(query_file=q, level='variant')

    # Run all required analysis for each query
    def run_analysis(self, query_file, level):
        cutoffs = None
        # Check if this query requires cutoffs to be defined by average values
        if query_file in Queries.AVERAGE_STDDEV:
            prequery = Queries.AVERAGE_STDDEV[query_file]
            cutoffs = self.average_stddev_cutoffs(prequery)
        query = self.prepare_query(query_file, preset_cutoffs=cutoffs)
        result = self.bq.run(query)
        failed = []
        if level == 'sample':
            failed = self.get_failed_samples(result)
        elif level == 'variant':
            failed = []
        return failed

    #### Query set up ####
    # Set up the query, read it in, apply substitutions
    def prepare_query(self, query_file, preset_cutoffs=None):
        raw_query = self.get_query(query_file)
        if preset_cutoffs is None:
            preset_cutoffs = self.get_preset_cutoffs(query_file)
        query = self.query_substitutions(raw_query, other=preset_cutoffs)
        return query

    # Read raw query in from file
    def get_query(self, file):
        path = os.path.join(self.query_repo, file)
        query = ''
        with open (path, "r") as f:
            query = f.read()
        return query

    # Apply any substitutions. Substitutions set on other must be in a dictionary
    def query_substitutions(self, query, other=None):
        replacements = {
            "_THE_TABLE_": Config.VARIANT_TABLE,
            "_THE_EXPANDED_TABLE_": Config.EXPANDED_TABLE,
            "_PATIENT_INFO_": Config.PATIENT_INFO,
            "_GENOTYPING_TABLE_": Config.GENOTYPING_TABLE
        }
        for r in replacements:
            query = query.replace(r, replacements[r])
        if other is not None:
            for r in other:
                query = query.replace(r, other[r])
        return query

    # Get preset cutoffs from query file
    def get_preset_cutoffs(self, query):
        cutoffs = None
        if query in Queries.PRESET_CUTOFFS:
            cutoffs = Queries.PRESET_CUTOFFS[query]
        return cutoffs

    # Run metrics query to define cutoffs based on average and standard deviation values
    def average_stddev_cutoffs(self, query_file):
        query = self.prepare_query(query_file)
        result = self.bq.run(query)
        average, stddev = self.get_average_stddev(result)
        max, min = self.calculate_max_min(average, stddev)
        substitutions = self.create_max_min_substitutions(max, min)
        return substitutions

    # Get average and standard deviation values from a parsed BigQuery result
    def get_average_stddev(self, result):
        for r in result:
            average = r['average']
            stddev = r['stddev']
            return average, stddev

    # Calculate maximum and minimum values based on average and standard deviation.
    # Cutoffs are defined as more than three standard deviations away from the average.
    def calculate_max_min(self, average, stddev):
        max = float(average) + (3 * float(stddev))
        min = float(average) - (3 * float(stddev))
        return max, min

    # Create maximum and minimum cutoff dictionary
    def create_max_min_substitutions(self, max, min):
        dict = {
            "_MAX_VALUE_": "%s" % max,
            "_MIN_VALUE_": "%s" % min,
        }
        return dict

    # Get failed sample ids
    def get_failed_samples(self, result):
        if result:
            failed_ids = []
            for r in result:
                failed_ids.append(r['sample_id'])
            return failed_ids
        return None

    # Get failed positions
    def get_failed_positions(self, result):
        if result:
            failed_positions = []
            for r in result:
                position = {
                    "chr": r['reference_name'],
                    "start": r['start'],
                    "end": r['end']
                }
                failed_positions.append(position)
            return failed_positions
        return None

    #### Miscellaneous functions ####
    # Check if a path exists
    def check_path(self, file):
        if not os.path.exists(file):
            raise Exception("%s not found!" % file)

    def setup_log(self, verbose):
        log_level = ''
        if verbose is True:
            log_level = logging.DEBUG
        else:
            log_level = logging.INFO
        # todo write to file
        logging.basicConfig(format='%(asctime)s %(message)s', level=log_level)


qc = GenomicsQC(verbose=True)
qc.variant_qc()
#qc.sample_qc()
