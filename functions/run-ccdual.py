# Run results: use python os module
# Date: Oct 9, 2011

import os
import time
instancepath = "../instances/"
resultpath = "../results/"

command = []

instance = ["1-4-multi"];
scenario = ["500"];
resultfile = "temp";
option =  [" 1 2 1 0.001 "]

for m in instance:
    for n in scenario:
        for k in option:
            for i in range(5):
                command = command + ["./ccdual-decomp " + instancepath + m + "-" + n + "-" + str(i+1) + ".txt " + " 0.1 " + k + resultpath + resultfile];
                command = command + ["./ccdual-decomp " + instancepath + m + "-" + n + "-" + str(i+1) + ".txt " + " 0.2 " + k + resultpath + resultfile];

for i in range(len(command)):
	print command[i];
	os.system(command[i]);
	
