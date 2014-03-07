#!/usr/bin/python

from subprocess import call
import time
import sys
from datetime import timedelta
import datetime
start = time.time()

eta=float(sys.argv[1])
eta_max=float(sys.argv[2])
step=float(sys.argv[3])

while (eta < eta_max):
	print ("Starting sfere2D with eta=" + str(eta))
	call(["./main",str(eta)])
	eta=eta+step

print("Script exited successfully")
print("time elapsed: ", str(datetime.timedelta(seconds=(time.time() - start))))