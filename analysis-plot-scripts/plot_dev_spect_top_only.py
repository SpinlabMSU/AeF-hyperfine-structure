#!/usr/bin/env python3
## plot_dev_spect.py -- makes plots with the same format as figure (3) from
## PRA 98, 032513 for devonshire-potential enabled 138BaF systems. 
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

import matplotlib


matplotlib.rcParams.update({'font.size': 20})

nocolor = False
plottype = 'internal'
nobox = False
relative_energy = False

color_pz = 'b'
color_xy = 'g'
color_nz = 'r'

#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-185338.4873197' #nodev
#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779' #deven
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\test\2024-04-05-232834.8866135' # good embedded-in-medium run
#r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-12-193933.6112353'

if len(sys.argv) > 1:
    rundir = sys.argv[1]

if len(sys.argv) > 2:
    plottype = sys.argv[2]


if plottype == 'mcaw24':
    nocolor = True
    nobox = True
    relative_energy = True
#
if nocolor:
    color_pz = 'black'
    color_xy = 'black'
    color_nz = 'black'

#
rundir = os.path.abspath(rundir)
run = os.path.split(rundir)[1]

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]
print(Ez)
print(states)
Ezs = df[Ez]
mid_idx = len(Ezs) // 2
Ez_mid = Ezs[mid_idx] / 1000

do_extras = False
# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

##
def get_energies(n):
    return df[f'E{n}']
## 
## suspect MDA is actually backwards
# +Z oriented states start at 0
bidx_pz = 0
Ep0s = get_energies(bidx_pz + 0)
Ep1s = get_energies(bidx_pz + 1)
Ep2s = get_energies(bidx_pz + 2)
Ep3s = get_energies(bidx_pz + 3)

Em0s = get_energies(bidx_pz + 4)
# -z oriented states
bidx_nz = 20 # 4 (for -z) + 4*4 (ffor +-x, +-y)
En0s = get_energies(bidx_nz + 0)
En1s = get_energies(bidx_nz + 1)
En2s = get_energies(bidx_nz + 2)
En3s = get_energies(bidx_nz + 3)

delta_10_mf0 = Ep1s[1] - Ep0s[1]
delta_10_mf1 = Ep2s[1] - Ep0s[1]
delta_10_mft = Ep3s[1] - Ep0s[1]

delta_10_pmf1 = Ep2s[1] - Ep0s[1]
delta_10_pmft = Ep3s[1] - Ep0s[1]
delta_10_nmf1 = En2s[1] - En0s[1]
delta_10_nmft = En3s[1] - En0s[1]

dn0 = Ep1s[1] - Ep0s[1]#np.average(En1s - En0s)
dp0 = En1s[1] - En0s[1]#np.average(Ep1s - Ep0s)

if relative_energy:
    Ep0s -= Em0s
    En0s -= Em0s
    Em0s -= Em0s

### Make equivalent to PRA fig 3.a
props = dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.5)
multiplier = 5
x_pix = 329 * multiplier
y_pix = 173 * multiplier
fontsize = 28
fig,ax = plt.subplots(1, 1, figsize=(x_pix / 100.0, y_pix / 100.0), sharex = True, sharey = False)
# part a
gca = plt.gca() #plt.subplot(4, 1, 1)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'Stark Shift of Lowest-energy group of states'
if not nobox:
    gca.text(0.015, 0.12, textstr, transform=gca.transAxes, fontsize=28, verticalalignment='top', bbox=props)


color = []
for i in range(4): color.append(color_pz)
for i in range(16): color.append(color_xy)
for i in range(4): color.append(color_nz)

Ezs_kV = np.array(df[Ez]) / 1000.0
#df.plot(Ez[:24], states[:24], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = False)
if not relative_energy:
    plt.ylabel('Energy (GHz)')
else:
    plt.ylabel('Relative Stark Shift (GHz)')
#plt.plot(df[Ez][1:], df[states[:24]][1:] / 1000)
# +Z
plt.plot(Ezs_kV, Ep0s/1000, color=color_pz)
plt.annotate('$+\hat{Z}$', xy=(Ez_mid, Ep0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 1.5), color=color_pz, textcoords='offset points')
# XY
plt.plot(Ezs_kV, Em0s/1000, color=color_xy)
plt.annotate('$+\hat{X},-\hat{X},+\hat{Y},-\hat{Y}$', xy=(Ez_mid, Em0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 2.5), color=color_xy, textcoords='offset points')
# -Z
plt.plot(Ezs_kV, En0s/1000, color=color_nz)
plt.annotate('$-\hat{Z}$', xy=(Ez_mid, En0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 11.5), color=color_nz, textcoords='offset points')
gca.set_xlim([0, 50])

#fig.suptitle(f"N=0, F=0,1 Stark shift for run {run}", y=0.999)
plt.subplots_adjust(bottom=0.09, right=0.990, top=0.999, left = 0.09, hspace = 0.0)
plt.xlabel("Externally-Applied Electric field Strength (kV/cm)")

plt.savefig(os.path.join(rundir, f'pra_aspect_deven_top_{plottype}.png'))
plt.show()
sys.exit(0)