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




files=[]

file_temp = glob.glob('checkedpression*.dat')	
for f in file_temp:
	os.rename(f,f[7:])
