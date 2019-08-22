import pandas
from googleapiclient import discovery, http
from google_auth_oauthlib import get_user_credentials
from os import path
from re import sub
from io import FileIO
from slugify import slugify
from sys import exit
import json
import logging
from email.mime.text import MIMEText
import base64



###############################
#### FUNCTIONS DEFINITIONS ####
###############################

def get_proper_mimetype_and_filename(Mod_time: str, name: str, mimetype: str):
    # Here we change names and such, so that they are of correct format for file name. slugify is amazing at this
    filenamed_name = slugify(name)
    Mod_time_no_wierd_symbols = sub(':|\.', '_', Mod_time)

    # Here we establish which files should be saved on disc using which format, based on their mimetype. This needs to be inside the function, cause it establishes two
    mimetype_to_export_file_as = 'text/plain'
    extension_for_exported_file = '.txt'

    if mimetype == 'application/vnd.google-apps.spreadsheet':
        mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        extension_for_exported_file = '.xlsx'
    elif mimetype == 'application/vnd.google-apps.document':
        mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
        extension_for_exported_file = '.docx'
    elif mimetype == 'application/vnd.google-apps.presentation':
        mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.presentationml.presentation'
        extension_for_exported_file = '.pptx'

    return {'proper_mimetype' : mimetype_to_export_file_as, 'proper_filename' : f'{filenamed_name}_{Mod_time_no_wierd_symbols}{extension_for_exported_file}'}



# The parameters are the data you get from google drive
def download_the_file(Id: str, Mod_time: str, name: str, mimetype: str, drive):

    proper_mimetype_and_filename = get_proper_mimetype_and_filename(Mod_time = Mod_time, name = name, mimetype = mimetype)

    # I am not sure why this workflow looks the way it looks. Future me - your on your own and you know as much as I do
    try:
        request = drive.files().export(fileId = Id, mimeType = proper_mimetype_and_filename['proper_mimetype'])
    except:
        logging.error('creating export failed')
        raise Exception('xxx')


    try:
        fh = FileIO(proper_mimetype_and_filename['proper_filename'], 'wb')
    except:
        logging.error('creating file on a disc failed')
        raise Exception('xxx')

    downloader = http.MediaIoBaseDownload(fh, request)

    done = False
    while done is False:
        try:
            status, done = downloader.next_chunk()
        except:
            logging.error('downloading piece of a file failed')
            raise Exception('xxx')

        print(f"Downloading: {proper_mimetype_and_filename['proper_filename']}")



# This function compares old andnew time
def download_the_file_with_changed_modification_time(Id: str, Mod_time_old: str, Mod_time_new: str, name: str, mimetype: str, drive):
    if Mod_time_old != Mod_time_new:
        download_the_file(Id = Id, Mod_time = Mod_time_new, name = name, mimetype = mimetype, drive = drive)



# The secret file needs to be in
def set_id_and_secret(secret_file_name_: str):
    if path.isfile(secret_file_name_) == True:
        with open(secret_file_name_) as json_file:
            json_file_loaded = json.load(json_file)
            print()
            return {'client_id' : json_file_loaded['installed']['client_id'], 'client_secret' : json_file_loaded['installed']['client_secret']}
    else:
        logging.error(f'No secret file found. It needs to be inside the same folder as the script and should be named: {secret_file_name_}')
        raise Exception('xxx')

#### WORK IN PROGRESS!! #### try send/import? Tutaj message ma złą składnie
def send_message(service, user_id, message):
    message = service.users().messages().send(userId=user_id, body=message).execute()
    return message

def create_message(sender: str = 'adrianstankiewicz85@gmail.com', to: str = 'adrianstankiewicz85@gmail.com', subject: str = 'asshole', message_text: str = 'asshole text'):
    message = MIMEText(message_text)
    message['to'] = to
    message['from'] = sender
    message['subject'] = subject

    # This is what need to be done to send mail via google api. Im human trash to them. MIMEText -> string -> bytes -> string(??). Cool.
    bytes_message = base64.urlsafe_b64encode(message.as_string().encode())
    bytes_message_decoded = bytes_message.decode()
    return {'raw': bytes_message_decoded}
#### WORK IN PROGRESS!! ####




###############################
#### FUNCTIONS DEFINITIONS ####
###############################



def main():

    logging.basicConfig(format = '%(asctime)s - %(message)s', filename = 'GDocs_backuper.log')

    scopes = ['https://www.googleapis.com/auth/drive', 'https://www.googleapis.com/auth/gmail.send']
    secret_file_name = 'client_secret.json'

    client_id = set_id_and_secret(secret_file_name_ = secret_file_name)['client_id']
    client_secret = set_id_and_secret(secret_file_name_ = secret_file_name)['client_secret']

    try:
        creds = get_user_credentials(scopes, client_id, client_secret)
    except:
        logging.error('get_user_credentials failed')
        raise Exception('get_user_credentials failed')


    # Construct a Resources/services for interacting with an API. Building a Python representation of the API
    try:
        drive = discovery.build(serviceName = 'drive', version = 'v3', credentials = creds)
        mail = discovery.build(serviceName = 'gmail', version = 'v1', credentials = creds)
    except:
        logging.error('building a drive or mail image failed')
        raise Exception('building a drive or mail image failed')

    #### WORK IN PROGRESS!! ####
    did_u_sent_it = send_message(service = mail, user_id = 'me', message = create_message())
    #### WORK IN PROGRESS!! ####

    # execute() sends the request to google. The things before it, are methods for this specific resource type. The "get" is a python dictionary method. Not sure how Resoucre objects translates to python dictionary though. Dont understand get() arguments either. list() needs to explicitly state needed fields, or it returns only small set of them. "*" returns all fields. You can return only fields You want, but syntax of this query was written by some cunt, and its a waste of my time to learn it.
    try:
        files = drive.files().list(fields = '*').execute().get('files', [])
    except:
        logging.error('getting list of files failed')
        raise Exception('xxx')


    # list of files of interest. This is wierd way to make dataframes in padnas work...
    new_files_of_interest = pandas.DataFrame(data = {'ID': [], 'ModifiedTime': [], 'Name': [], 'MimeType': []}, dtype = str)

    # Here we write a dataframe of google documents in the drive currently
    for file in files:
        if file['mimeType'] in ['application/vnd.google-apps.document', "application/vnd.google-apps.spreadsheet", "application/vnd.google-apps.presentation"]: # So for now lets return only google docs. I dunno if we want to backup everything there is on the drive
            new_files_of_interest = new_files_of_interest.append( other = pandas.DataFrame(data = {'ID': [file['id']], 'ModifiedTime': [file['modifiedTime']], 'Name': [file['name']], 'MimeType': [file['mimeType']]}, dtype = str) ) #



    # Here we either download all google documents, if this is a first time the script is running, or we identify which documents chages since last time, and we download only those. It would be nice to limit number of documents to 100 or 1000 and sign them with data of download
    if path.isfile('files_of_interest.csv') == True:
        old_files_of_interest = pandas.read_csv('files_of_interest.csv')

        combined_files_of_interest = new_files_of_interest.merge(right = old_files_of_interest, how = 'left', on = 'ID', suffixes = ('_new', '_old'))

        try:
            new_files_of_interest.to_csv('files_of_interest.csv')
        except:
            logging.error('creating new file list failed')
            raise Exception('xxx')

        ### Run this to finalize the file download
        combined_files_of_interest.apply(
            func = lambda a : download_the_file_with_changed_modification_time(Id = a.ID, Mod_time_old = a.ModifiedTime_old, Mod_time_new = a.ModifiedTime_new, name = a.Name_new, mimetype = a.MimeType_new, drive = drive),
            axis = 1)
    else:
        # Here we save the first list of google doc files on the drive
        try:
            new_files_of_interest.to_csv('files_of_interest.csv')
        except:
            logging.error('creating new file list failed')
            raise Exception('xxx')


        # Here we download all files in file list
        new_files_of_interest.apply(
            func = lambda a : download_the_file(Id = a.ID, Mod_time = a.ModifiedTime, name = a.Name, mimetype = a.MimeType, drive = drive),
            axis = 1)



if __name__ == "__main__":
	main()
