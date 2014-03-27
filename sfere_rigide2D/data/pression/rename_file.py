#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import sys
import glob
from operator import itemgetter
import os
import math



N = int(sys.argv[1])

files=[]

file_temp = glob.glob('%s/checkedpression*.dat'%(N))	
for f in file_temp:
	os.rename(f,f[:4]+f[11:])
