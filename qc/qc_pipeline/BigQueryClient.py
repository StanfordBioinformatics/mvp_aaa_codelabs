import httplib2
from apiclient.discovery import build
from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client import tools

class BigQueryClient(object):
    def __init__(self, client_secrets):
        self.client_secret_path = client_secrets

    def bigquery_setup(self):
        FLOW = flow_from_clientsecrets(self.client_secret_path,
                                   scope='https://www.googleapis.com/auth/bigquery')

        storage = Storage('bigquery_credentials.dat')
        credentials = storage.get()

        if credentials is None or credentials.invalid:
            credentials = tools.run_flow(FLOW, storage, tools.argparser.parse_args([]))

        http = httplib2.Http()
        http = credentials.authorize(http)

        bigquery_service = build('bigquery', 'v2', http=http)
        return bigquery_service

class BigQuery(object):
    def __init__(self, project_number, client_secrets):
        client = BigQueryClient(client_secrets)
        self.service = client.bigquery_setup()
        self.project_number = project_number

    def run(self, query):
        response = self.execute_query(query)
        if self.check_response(response) is False:
            return None
        result = self.parse_bq_response(response)
        return result

    # Run a query on BigQuery
    def execute_query(self, query):
        query_dict = {'query': query, 'timeoutMs': 600000}
        query_request = self.service.jobs()
        query_response = query_request.query(projectId=self.project_number, body=query_dict).execute()
        return query_response

    # Parse BigQuery response into list of dictionaries
    def parse_bq_response(self, response):
        fields = self.get_fields(response)
        result = []
        if 'rows' in response:
            for row in response['rows']:
                result_row = {}
                position = 0
                for field in row['f']:
                    result_row[fields[position]] = field['v']
                    position += 1
                result.append(result_row)
            return result
        return None

    # Get field names from BigQuery response
    def get_fields(self, response):
        fields = {}
        if 'schema' in response:
            position = 0
            for field in response['schema']['fields']:
                fields[position] = field['name']
                position += 1
            return fields
        return None

    # Check for valid response
    # todo add more checks
    def check_response(self, response):
        if 'jobComplete' in response:
            if response['jobComplete'] is True:
                return True
        raise False