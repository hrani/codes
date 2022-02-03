#listing all the python files in the directory and running those files
import os
from subprocess import call
directory = '/home/rani/moose-examples/tutorials/'
listOfFiles = list()
for (dirpath, dirnames, filenames) in os.walk(directory):
    listOfFiles += [os.path.join(dirpath, file) for file in filenames]
#print(" list ", listOfFiles)
for filename in listOfFiles:
    if filename.endswith(".py"):
        print("FN ",filename)
        if filename not in ["multiscaleOneCompt.py"]:
            ts =  os.path.join(directory, filename)
            call(["python3", ts])
            continue
    else:
        continue
'''
for filename in os.listdir(directory):
    if filename.endswith(".py"):
        print("FN ",filename)
        if filename not in ["multiscaleOneCompt.py"]:
            ts =  os.path.join(directory, filename)
            call(["python3", ts])
            continue
    else:
        continue
'''
