from __future__ import print_function

from googleapiclient import discovery
from httplib2 import Http
import google_auth_oauthlib
from google_auth_oauthlib.flow import Flow, InstalledAppFlow
from oauth2client import file, client, tools

scopes = 'https://www.googleapis.com/auth/drive.metadata.readonly'
client_id = '547846840425-0f63erbq9s1vtr9htm6j5l8d4ujr8gck.apps.googleusercontent.com'
client_secret = 'tMCbjTRP5FStZCOEfS2EMKj1'
# store = file.Storage('storage.json')

# flow = InstalledAppFlow.from_client_secrets_file(client_secrets_file = 'client_secret.json',  scopes = SCOPES)



# # Create the flow using the client secrets file from the Google API Console. Just an instance of a class, I think.
# flow = Flow.from_client_secrets_file(client_secrets_file = 'client_secret.json', scopes = scopes, redirect_uri = 'urn:ietf:wg:oauth:2.0:oob')
#
# # Generates an authorization URL. Tell the user to go to the authorization URL. I dont know what "prompt: consent" means, but it was inside authorization_url method
# niggers, _ = flow.authorization_url()
# print('Go here, you fucking whore: {}'.format(niggers))
#
# # The user will get an authorization code. This code is used to get the access token.
# code = input('Enter the authorization code: ')
# flow.fetch_token(code = code)
#
# # You can use flow.credentials, or you can just get a requests session using flow.authorized_session.
# # session = flow.authorized_session()
# creds = flow.credentials
creds = google_auth_oauthlib.get_user_credentials(scopes, client_id, client_secret)


# Construct a Resource for interacting with an API.
DRIVE = discovery.build(serviceName = 'drive', version = 'v3', credentials = creds)

files = DRIVE.files().list().execute().get('files', [])
for f in files:
    print(f['name'], f['mimeType'])
