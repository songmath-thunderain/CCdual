# Run results: use python os module
# Date: Oct 9, 2011

import os
import time
instancepath = "../instances/"
resultpath = "../results/"

command = []

instance = ["sc1-100"];
scenario = ["1000", "5000", "10000"];
resultfile = "Dec14-2015-coveringLP";
option =  [" 4 2 2 0.001 "]

for m in instance:
    for n in scenario:
        for k in option:
            for i in range(5):
                command = command + ["./ccdual-decomp2 " + instancepath + m + "-" + n + "-" + str(i+1) + ".txt " + " 0.1 " + k + resultpath + resultfile + " 0 "];

for i in range(len(command)):
	print command[i];
	os.system(command[i]);
	
