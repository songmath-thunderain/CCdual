# Run results: use python os module
# Date: Oct 9, 2011

import os
import time
instancepath = "../instances/"
resultpath = "../results/"

command = []

instance = ["1-4-multi", "1-6-multi"];
scenario = ["100", "500", "1000"];
resultfile = "temp";

for m in instance:
    for n in scenario:
            for i in range(5):
                command = command + ["./ccdual-bounds " + instancepath + m + "-" + n + "-" + str(i+1) + ".txt " + " 0.1 "  + resultpath + resultfile];

for m in instance:
    for n in scenario:
            for i in range(5):
                command = command + ["./ccdual-bounds " + instancepath + m + "-" + n + "-" + str(i+1) + ".txt " + " 0.2 "  + resultpath + resultfile];

for i in range(len(command)):
	print command[i];
	os.system(command[i]);
	
