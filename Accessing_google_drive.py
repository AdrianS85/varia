from __future__ import print_function

import pandas
from googleapiclient import discovery, http
from google_auth_oauthlib import get_user_credentials
from os import path
from re import sub
from io import FileIO
from slugify import slugify



def main():
    ###############################
    #### FUNCTIONS DEFINITIONS ####
    ###############################

    # The parameters are the data you get from google drive
    def Download_the_file(Id: str, Mod_time: str, name: str, mimetype: str):
        # Here we change names and such, so that they are of correct format for file name. slugify is amazing at this
        name_no_spaces = slugify(name)
        Mod_time_no_wierd_symbols = sub(':|\.', '_', Mod_time)

        # Here we establish which files should be saved on disc using which format, based on their mimetype
        mimetype_to_export_file_as = None
        if mimetype == 'application/vnd.google-apps.spreadsheet':
            mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        elif mimetype == 'application/vnd.google-apps.document':
            mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
        elif mimetype == 'application/vnd.google-apps.presentation':
            mimetype_to_export_file_as = 'application/vnd.openxmlformats-officedocument.presentationml.presentation'
        else:
            mimetype_to_export_file_as = 'text/plain'

        # I am not sure why this workflow looks the way it looks. Future me - your on your own and you know as much as I do
        request = DRIVE.files().export(fileId = Id, mimeType = mimetype_to_export_file_as)
        print('Tring to download: ', f"Downloading: f'{name_no_spaces}_{Mod_time_no_wierd_symbols}.goo'")
        fh = FileIO(f'{name_no_spaces}_{Mod_time_no_wierd_symbols}.goo', 'wb')
        downloader = http.MediaIoBaseDownload(fh, request)
        done = False
        while done is False:
            status, done = downloader.next_chunk()
            print(f"Downloading: f'{name_no_spaces}_{Mod_time_no_wierd_symbols}.goo'")



    # This function compares old andnew time
    def Compare_Times(Id: str, Mod_time_old: str, Mod_time_new: str, name: str, mimetype: str):
        if Mod_time_old != Mod_time_new:
            Download_the_file(Id = Id, Mod_time = Mod_time_new, name = name, mimetype = mimetype)
        else:
            pass

    ###############################
    #### FUNCTIONS DEFINITIONS ####
    ###############################



    ###############################
    ########## SCRIPT #############
    ###############################

    scopes = 'https://www.googleapis.com/auth/drive'
    client_id = 'xxx'
    client_secret = 'xxx'

    creds = get_user_credentials(scopes, client_id, client_secret)

    # Construct a Resource for interacting with an API. Building a Python representation of the API
    DRIVE = discovery.build(serviceName = 'drive', version = 'v3', credentials = creds)

    # execute() sends the request to google. The things before it, are methods for this specific resource type. The "get" is a python dictionary method. Not sure how Resoucre objects translates to python dictionary though. Dont understand get() arguments either. list() needs to explicitly state needed fields, or it returns only small set of them. "*" returns all fields. You can return only fields You want, but syntax of this query was written by some cunt, and its a waste of my time to learn it.
    files = DRIVE.files().list(fields = '*').execute().get('files', [])



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
        # combined_files_of_interest.to_csv('files_of_interest.csv')

        ### Run this to finalize the file download
        combined_files_of_interest.apply(
            func = lambda a : Compare_Times(Id = a.ID, Mod_time_old = a.ModifiedTime_old, Mod_time_new = a.ModifiedTime_new, name = a.Name_new, mimetype = a.MimeType_new),
            axis = 1)
    else:
        # Here we save the first list of google doc files on the drive
        new_files_of_interest.to_csv('files_of_interest.csv')

        # Here we download all files in file list
        new_files_of_interest.apply(
            func = lambda a : Download_the_file(Id = a.ID, Mod_time = a.ModifiedTime, name = a.Name, mimetype = a.MimeType),
            axis = 1)

    ###############################
    ########## SCRIPT #############
    ###############################

if __name__ == "__main__":
	main()
