#!/usr/bin/env python3
## plot_dev_spect.py -- makes plots with the same format as figure (2) from
## PRA 98, 032513 for devonshire-potential disabled (vacuum) 138BaF systems. 
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
import scipy.stats as scistat
import scipy.optimize as spopt
import aef_run


rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779'
if len(sys.argv) > 1:
    rundir = sys.argv[1]

rundir = os.path.abspath(rundir)

run = aef_run.aef_run(rundir)
#plt.rcParams['text.usetex'] = True
plt.rcParams["mathtext.fontset"] = 'cm'

x_pix = 1366
y_pix = 900
doe_plot = True

if len(sys.argv) > 3:
    x_pix = int(sys.argv[2])
    y_pix = int(sys.argv[3])

out_fname = "stark_shift_gnd.png"

# default fname for pra aspect ratios
if x_pix == 756:# and y_pix == 857:
    out_fname = "stark_shift_gnd_pra_aspect.png"

if len(sys.argv) > 4:
    out_fname = sys.argv[4]
run = os.path.split(rundir)[1]

starkpath = os.path.join(rundir, 'stark_shift_gnd.csv')
df = pd.read_csv(starkpath)

## csv keys -- note that the spaces in front are intentional
key_E = 'E-field (V/cm)'
key_dEgnd = ' dE_gnd'
key_dE0 = ' dE_0'
key_dE1 = ' dE_1'
key_dE2 = ' dE_2'
key_dE3 = ' dE_3'
key_dE4 = ' dE_4'
key_dE5 = ' dE_5'
key_dE6 = ' dE_6'
key_dE7 = ' dE_7'


df[key_E] /= 1000
Ezs = df[key_E]
dE_gnds = df[key_dEgnd]
dE_f1ts = df[key_dE2]
dE_f10s = df[key_dE3]
dE_f11s = df[key_dE4]
mid_idx = len(Ezs) // 2
Ez_mid = Ezs[mid_idx]


## make linear fit
Ez_uhalf = Ezs[mid_idx:]
Eg_uhalf = dE_gnds[mid_idx:]
lresult = scistat.linregress(Ez_uhalf, Eg_uhalf)

def line(Ez):
    m = lresult.slope
    b = lresult.intercept
    return m * Ez + b

## make quadratic fit


Ezs_line = Ezs[mid_idx // 2:]
dE_gnd_linears = line(Ezs_line)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig,ax = plt.subplots(2, 1, figsize=(x_pix / 100.0, y_pix / 100.0), sharex = True, sharey = False)
## Overall shift
gca = plt.subplot(2, 1, 1)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'Ground State Stark Shift'
gca.text(0.05, 0.1, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
plt.plot(Ezs, dE_gnds / 1000, 'b-')
plt.plot(Ezs_line, dE_gnd_linears / 1000, 'r--')
plt.ylabel("Stark Shift (GHz)")
#df.plot(key_E, [key_dEgnd,], ax=gca, ylabel='Energy (MHz)')
# Ground-relative shift
gca = plt.subplot(2, 1, 2)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'$F_1=1$ Stark shift above ground'
gca.text(0.05, 0.1, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
df.plot(key_E, [key_dE2, key_dE3, key_dE4], ylabel='Energy above $f_1$=0 ground states (MHz)', ax=gca, legend=False)
plt.xlabel("Externally-applied electric field strength (kV/cm)")
plt.annotate('$f=\\frac{3}{2}$', xy=(Ez_mid, dE_f11s[mid_idx - 1]), xycoords='data', xytext=(1.5, -12.5), color='g', textcoords='offset points', fontsize = 14)
plt.annotate('$f=\\frac{1}{2}$', xy=(Ez_mid, dE_f10s[mid_idx - 1]), xycoords='data', xytext=(1.5, 1.5), color='orange', textcoords='offset points', fontsize = 14)
#
if not doe_plot:
    fig.suptitle(f"N=0, F$_1$=0,1 Stark shift for run {run}", y=0.999)
else:
    fig.suptitle(f"N=0, F$_1$=0,1 Stark shift of Gas-phase heavy-half with no hyperfine structure on the light nucleus", y=0.999)
plt.subplots_adjust(bottom=0.05, right=0.95, top=0.97, left = 0.114, hspace = 0.0)
# 
plt.savefig(os.path.join(rundir, out_fname), bbox_inches="tight")
plt.show()