#listing all the python files in the directory and running those files
import os
from subprocess import call
directory = '/home/harsharani/moose/hrani_moose-example_cleanuppatch/snippets/'
for filename in os.listdir(directory):
    if filename.endswith(".py"):
        ts =  os.path.join(directory, filename)
        call(["python3", ts])
        continue
    else:
        continue
