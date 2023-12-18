#!/usr/bin/env python3
## plot_spectrum.py -- plots the spectrum of the hamiltonian as a function of
## the externally-applied electric field. 
# This file is part of the AeF-hyperfine-structure program. 
    
# AeF-hyperfine-structure is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version.

# AeF-hyperfine-structure is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

# You should have received a copy of the GNU General Public License along with
# AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
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
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779' #deven

if len(sys.argv) > 1:
    rundir = sys.argv[1]
run = os.path.split(rundir)[1]

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]

# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Energy Spectrum for run {run}")
df.plot(Ez, states, ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
plt.show()
