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
            #Run oauth2 flow with default arguments.
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

    def execute_query(self, query):
        query_dict = {'query': query, 'timeoutMs': 600000}
        query_request = self.service.jobs()
        query_response = query_request.query(projectId=self.project_number, body=query_dict).execute()
        return query_response

#client_secret_path = '/Users/gmcinnes/bin/google_genomics/client_secrets.json'
#query = 'SELECT TOP(title, 10) as title, COUNT(*) as revision_count FROM [publicdata:samples.wikipedia] WHERE wp_namespace = 0;'
#PROJECT_NUMBER = '963911152157'
#print query
#bq = BigQuery(project_number=PROJECT_NUMBER, client_secrets=client_secret_path)
#response = bq.execute_query(query)
#print response
#for row in response['rows']:
#  result_row = []
#  for field in row['f']:
#    result_row.append(field['v'])
#  print ('\t').join(result_row)
#
#print response.rows
