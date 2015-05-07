import argparse
import httplib2
import pprint
from apiclient.discovery import build
from collections import Counter
from oauth2client import tools
from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client.tools import run_flow

class GoogleGenomicsClient(object):
    def __init__(self, client_secrets, project_number):
        self.client_secrets_path = client_secrets
        self.project_number = project_number
        self.service = self.setup(self.client_secrets_path)


    def setup(self, client_secrets):
        storage = Storage('credentials.dat')
        credentials = storage.get()
        if credentials is None or credentials.invalid:
            flow = flow_from_clientsecrets(
                client_secrets,
                scope='https://www.googleapis.com/auth/genomics',
                message='You need to copy a client_secrets.json file into this directory, '
                'or pass in the --client_secrets_filename option to specify where '
                'one exists. See the README for more help.')
            credentials = run_flow(flow, storage)
            # Create a genomics API service
        http = httplib2.Http()
        http = credentials.authorize(http)
        service = build('genomics', 'v1beta2', http=http)
        return service

    def datasets(self):
        request = self.service.datasets().list(projectNumber=self.project_number)
        response = request.execute()
        pprint.pprint(response)

#client_secret_path = '/Users/gmcinnes/bin/google_genomics/client_secrets.json'
#g = GoogleGenomicsClient(client_secrets=client_secret_path)
#g.datasets()