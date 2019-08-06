#!/usr/bin/env python3

import sys
import shutil
import csv
import re
# python3 sciNotRem.py3 $(ls) <- this is how to turn ls to arguments



### This part shows user which files are up for sciNotRem
print("##############################################################\nYou have uploaded these files for removing scientific notation:\n##############################################################\n")

# Here we make list of only all files inputed to scrip as arguement except this file, and files already done by this script
files_to_do = []
for a in sys.argv:
    if (a != sys.argv[0] and a.find("_SciNotRepl.") == -1):
        files_to_do.append(a)

# We print it as a separate step to make sure that the list was created correctly
for b in files_to_do:
    print(b)

sys.stderr.write("\n##############################################################\nCheck if these files are the ones You want to change!!\n##############################################################\n")
sys.stderr.flush()
### This part shows user which files are up for sciNotRem



### This parts makes user be sure, that files are correct
while True:
        s = input('Do You wish to proceed? (y/n)\n')
        if s == "n":
            print("OK. Exiting now.")
            sys.exit()
        elif s == "y":
            print("Good.\n")
            break
        else:
            print('Sorry, I only understand "y" and "n" answers. Try again.')
            continue

### This parts makes user be sure, that files are correct



for c in files_to_do:
    print("File is being updated: " + c)
    with open(c) as csvfile:
        # Create output file for writting updated data
        new_name_for_updated_file = c.replace(".", "_SciNotRepl.")
        pre_writeit = open(new_name_for_updated_file, "w")
        writeit = csv.writer(pre_writeit, delimiter=',')

        # Update data
        readit = csv.reader(csvfile, delimiter=',')
        for row in readit:
            updated_row = [w.replace( 'Start', str(float(10)) ) for w in row]
            writeit.writerow(updated_row)
        pre_writeit.close()
        print("Updated file written: " + new_name_for_updated_file + "\n")
