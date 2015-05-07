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
            cutoff = None
            if q in Queries.PRESET_CUTOFFS:
                cutoff = Queries.PRESET_CUTOFFS[q]
            query = self.get_query(q)
            response = self.run_query(query)
            self.parse_bq_response(response)

    def variant_qc(self):
        print "Running Variant Level QC"


    def get_query(self, file):
        path = os.path.join(self.query_repo, file)
        query = self.read_query(path)
        preset_cutoffs = self.get_preset_cutoffs(file)
        query = self.query_substitutions(query, other=preset_cutoffs)
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
        bq = BigQueryClient.BigQuery(client_secrets=self.client_secrets_path, project_number=self.project_number)
        response = bq.execute_query(query)
        return response
        # todo make sure queries finish

    def parse_bq_response(self, response):
        print response
        self.check_response(response)
        failed_ids = []
        if 'rows' in response:
            for row in response['rows']:
                result_row = []
                for field in row['f']:
                    result_row.append(field['v'])
                failed_ids.append(result_row[0])
            return failed_ids

    def check_response(self, response):
        if 'jobComplete' in response:
            if response['jobComplete'] is True:
                return True
        raise Exception("Query failed to complete")
        return False

qc = GenomicsQC()
qc.sample_qc()
