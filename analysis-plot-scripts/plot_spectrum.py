# system imports
import math
import cmath
import sys
import os
import os.path
import re
#
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npla
import pandas as pd
import numba


#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-185338.4873197' #nodev
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-190613.6465128' #deven
run = os.path.split(rundir)[1]
if len(sys.argv) > 1:
    rundir = sys.argv[1]

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]

fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Energy Spectrum for run {run}")
df.plot(Ez, states, ylabel = 'Energy (MHz)', ax=plt.gca())
plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
plt.show()
