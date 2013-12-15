from __future__ import with_statement
# numpy
import numpy as np
# matplotlib plotting module
import matplotlib.pyplot as plt
# matplotlib colormap module
import matplotlib.cm as cm
# needed for formatting Y axis
from matplotlib.ticker import FuncFormatter
# Matplotlib font manager
import matplotlib.font_manager as font_manager

with open('population.csv') as f

y = np.loadtxt(f)
