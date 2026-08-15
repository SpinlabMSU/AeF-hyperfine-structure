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
import aef_run


#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-185338.4873197' #nodev
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779' #deven

if len(sys.argv) > 1:
    rundir = sys.argv[1]
rundir = os.path.abspath(rundir)
run = aef_run.aef_run(rundir)
run_name = run.run

#starkpath = os.path.join(rundir, 'stark_spectrum.csv')
#df = pd.read_csv(starkpath)
df = run.parse_stark_spect()
field_type = 'electric' if not run.is_zeeman else 'magnetic'

Ez = df.keys()[1]
print(Ez)
states = df.keys()[2:]

zero_ground = False
zero_idx = 0
measure_deviation = False
use_volts = False
do_cut = False
do_hard_cut = False
black_dots = False
elabel = 'Energy'

ymax = None
ymin = None
max_idx_soft = None
max_idx_hard = None
scale = 'm'
# map to scale from Megahertz to 
scale_map = {
    'k' : (1E+03, 'kHz'),
    'm' : (1E+00, 'MHz'),
    'g' : (1E-03, 'GHz'),
    't' : (1E-06, 'THz')
}

plt.rcParams['font.size'] = 14
title = None
outname = None
plot_vline = False
vline_x = 0.0

## Actually parse arguments
for idx in range(2, len(sys.argv)):
    arg = sys.argv[idx]
    lrg = arg.lower()
    if arg.startswith('-z'): zero_ground = True
    if arg.startswith('-Z'): # zero out an arbitrary state
        zero_ground = True
        zero_idx = int(sys.argv[idx + 1])
        idx += 1 # skip next argument
    if lrg.startswith('-m'): measure_deviation = True
    if lrg.startswith('-v'): use_volts = True
    if arg.startswith('-c'): # "soft cut": sets window, lowercase c specifically
        do_cut = True
        max_idx_soft = int(sys.argv[idx + 1])
        idx += 1 # skip next argument
    if arg.startswith('-C'): #  "Hard cut": doesn't plot higher states, uppercase C specifically
        do_hard_cut = True
        max_idx_hard = int(sys.argv[idx + 1])
        idx += 1 # skip next argument
    if lrg.startswith('-s'):
        # scale -- TODO really implement
        scale = sys.argv[idx + 1]
        idx += 1
    if lrg.startswith('-b'): black_dots = True
    if lrg.startswith('-t'):
        title = sys.argv[idx + 1]
        idx += 1
    if lrg.startswith('-o'):
        outname = sys.argv[idx + 1]
        idx += 1
    if lrg.startswith('-vline'):
        plot_vline = True
        vline_x = float(sys.argv[idx + 1])
        idx += 1
# set default output filename:
if outname == None:
    outname = "spectrum_plot.png" if not black_dots else "spectrum_plot_no_state.png"
# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

# set scale factor and label
scale_factor, scale_label = scale_map[scale]
elabel = f'"Absolute" Energy ({scale_label})'

if do_hard_cut:
    states = states[:max_idx_hard]
    #states = []

if measure_deviation:
    print(df[states])
    zero_field_Es = np.array(df.iloc[0, 2:])
    print(zero_field_Es)
    for state in states:
        df[state] = df[state] - zero_field_Es
        print(df[state])

if zero_ground:
    print(df[states])
    zero_field_Es = df[f'E{zero_idx}']
    elabel = f'Excitation Energy ({scale_label}) above state #{zero_idx}'
    print(zero_field_Es)
    for state in states:
        df[state] = df[state] - zero_field_Es
        print(df[state])

# process y-axis scaling
df[states] *= scale_factor

# perform cut
if do_cut:
    #states = states[:maxn]
    state = states[max_idx_soft]
    print(f"Cutting plot at state #{max_idx_soft} = {state}")
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
if run.is_zeeman:
    df[Ez] *= 10000 # Tesla to Gauss
    df[Ez] *= 1000  # Gauss to milliGauss
    xlab = "Externally-applied magnetic field strength (mG)"
elif not use_volts: # only for electric field runs
    df[Ez] /= 1000 
    xlab = "Externally-applied electric field strength (kV/cm)"

plot_kind = "line" if not black_dots else "scatter"
color = None if not black_dots else "black"
# do plot
fig = plt.figure(figsize=(13.66, 9.00))
if title == None:
    title = f"Energy Spectrum for run {run_name}"
plt.title(title)
if not black_dots:
    df.plot(Ez, states, xlabel=xlab, ylabel = f'{elabel}', ax=plt.gca(), legend = use_legend)
else:
    df.plot(Ez, states, xlabel=xlab, ylabel = f'{elabel}', ax=plt.gca(), legend = use_legend, color=color, linestyle='', marker='o')
if ymax != None:
    plt.ylim(bottom=ymin, top=ymax)

if plot_vline:
    plt.axvline(vline_x, color='r')

plt.savefig(os.path.join(rundir, outname))
#if not black_dots:
#    plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
#else:
#    plt.savefig(os.path.join(rundir, 'spectrum_plot_no_state.png'))
plt.show()
