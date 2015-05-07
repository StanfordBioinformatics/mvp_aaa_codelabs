import os
import sys
import re
from GenomicsQueries import Queries
import BigQueryClient
from config import Config
from GoogleGenomicsClient import GoogleGenomicsClient

class GenomicsQC(object):
    def __init__(self):
        self.query_repo = Config.QUERY_REPO
        self.variant_table = Config.VARIANT_TABLE
        self.expanded_table = Config.EXPANDED_TABLE
        self.client_secrets_path = Config.CLIENT_SECRETS
        self.project_number = Config.PROJECT_NUMBER

    def sample_qc(self):
        print "Running Sample Level QC"
        sample_queries = Queries.SAMPLE_LEVEL_QC_QUERIES
        for q in sample_queries:
            print q
            self.run_analysis(q)

    def variant_qc(self):
        print "Running Variant Level QC"

    def run_analysis(self, query_file):
        # Check for prerequisites
        cutoffs = None
        if query_file in Queries.AVERAGE_STDDEV:
            prequery = Queries.AVERAGE_STDDEV[query_file]
            cutoffs = self.average_stddev_cutoffs(prequery)
        query = self.prepare_query(query_file, preset_cutoffs=cutoffs)
        response = self.run_query(query)
        result = self.parse_bq_response(response)
        failed = self.get_failed_ids(result)
        print failed

    def prepare_query(self, query_file, preset_cutoffs=None):
        raw_query = self.get_query(query_file)
        if preset_cutoffs is None:
            preset_cutoffs = self.get_preset_cutoffs(query_file)
        query = self.query_substitutions(raw_query, other=preset_cutoffs)
        return query

    def get_query(self, file):
        path = os.path.join(self.query_repo, file)
        query = self.read_query(path)
        return query

    def read_query(self, file):
        with open (file, "r") as f:
            query = f.read()
        return query

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

    def get_preset_cutoffs(self, query):
        cutoffs = None
        if query in Queries.PRESET_CUTOFFS:
            cutoffs = Queries.PRESET_CUTOFFS[query]
        return cutoffs

    def check_path(self, file):
        if not os.path.exists(file):
            raise Exception("%s not found!" % file)

    def run_query(self, query):
        #print query
        bq = BigQueryClient.BigQuery(client_secrets=self.client_secrets_path, project_number=self.project_number)
        response = bq.execute_query(query)
        return response
        # todo make sure queries finish

    def parse_bq_response(self, response):
        print response
        self.check_response(response)
        result = []
        if 'rows' in response:
            for row in response['rows']:
                result_row = []
                for field in row['f']:
                    result_row.append(field['v'])
                result.append(result_row)
            return result
        return None

    def get_failed_ids(self, result):
        if result is None:
            return None
        failed_ids = []
        for r in result:
            failed_ids.append(r[0])
        return failed_ids

    def check_response(self, response):
        if 'jobComplete' in response:
            if response['jobComplete'] is True:
                return True
        raise Exception("Query failed to complete")
        return False

    def average_stddev_cutoffs(self, query_file):
        query = self.prepare_query(query_file)
        response = self.run_query(query)
        result = self.parse_bq_response(response)
        average, stddev = self.get_average_stddev(result)
        max, min = self.calculate_max_min(average, stddev)
        substitutions = self.create_max_min_substitutions(max, min)
        return substitutions

    def get_average_stddev(self, result):
        for r in result:
            average = r[0]
            stddev = r[1]
            return average, stddev

    def calculate_max_min(self, average, stddev):
        max = float(average) + (3 * float(stddev))
        min = float(average) - (3 * float(stddev))
        return max, min

    def create_max_min_substitutions(self, max, min):
        dict = {
            "_MAX_VALUE_": "%s" % max,
            "_MIN_VALUE_": "%s" % min,
        }
        return dict

qc = GenomicsQC()
qc.sample_qc()
