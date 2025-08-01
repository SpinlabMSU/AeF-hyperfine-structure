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
rundir = os.path.abspath(rundir)
run = os.path.split(rundir)[1]

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]

zero_ground = False
measure_deviation = False
use_volts = False
do_cut = False

ymax = None
ymin = None
max_idx = None
scale = 'm'
# map to scale from Megahertz to 
scale_map = {
    'k' : (1E+03, 'kHz'),
    'm' : (1E+00, 'MHz'),
    'g' : (1E-03, 'GHz'),
    't' : (1E-06, 'THz')
}


## Actually parse arguments
for idx in range(2, len(sys.argv)):
    arg = sys.argv[idx]
    lrg = arg.lower()
    if lrg.startswith('-z'): zero_ground = True
    if lrg.startswith('-m'): measure_deviation = True
    if lrg.startswith('-v'): use_volts = True
    if lrg.startswith('-c'):
        do_cut = True
        max_idx = int(sys.argv[idx + 1])
        idx += 1 # skip next argument
    if lrg.startswith('-s'):
        # scale -- TODO really implement
        scale = sys.argv[idx + 1]
        idx += 1

# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

if measure_deviation:
    print(df[states])
    zero_field_Es = np.array(df.iloc[0, 2:])
    print(zero_field_Es)
    for state in states:
        df[state] = df[state] - zero_field_Es
        print(df[state])

if zero_ground:
    print(df[states])
    zero_field_Es = df['E0']
    print(zero_field_Es)
    for state in states:
        df[state] = df[state] - zero_field_Es
        print(df[state])

# process y-axis scaling
scale_factor, scale_label = scale_map[scale]
df[states] *= scale_factor

# perform cut
if do_cut:
    #states = states[:maxn]
    state = states[max_idx]
    print(f"Cutting plot at state #{max_idx} = {state}")
    Es = np.array(df[state])
    Egs = np.array(df['E0'])
    ymax = np.max(Es)
    ymin = np.min(Egs)
    dy = ymax - ymin
    pct = 0.01
    ymax += pct * dy
    ymin -= pct * dy
    print(Es)
    print(f"ymin is {ymin}, ymax is {ymax}")

xlab = "Externally-applied electric field strength (V/cm)"
if not use_volts:
    df[Ez] /= 1000 
    xlab = "Externally-applied electric field strength (kV/cm)"

# do plot
fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Energy Spectrum for run {run}")
df.plot(Ez, states, xlabel=xlab, ylabel = f'Energy ({scale_label})', ax=plt.gca(), legend = use_legend)
if ymax != None:
    plt.ylim(bottom=ymin, top=ymax)
plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
plt.show()
