#!/usr/bin/env python3

import sys
import shutil
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
            print("Good.")
            break
        else:
            print('Sorry, I only understand "y" and "n" answers. Try again.')
            continue

### This parts makes user be sure, that files are correct

copied_files_to_do = []
for a in files_to_do:
    shutil.copyfile(a, a.replace(".", "_SciNotRepl."))
    copied_files_to_do.append(a.replace(".", "_SciNotRepl."))

for c in copied_files_to_do:
    print("New file copied: " + c)
    file_object = open(c, "r")

    file_object.close()


#


# DO NOT CHANGE THE FILES IN PLACE!! What if You give some other files as an input - it will kill them
