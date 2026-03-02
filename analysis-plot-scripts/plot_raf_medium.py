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
import aef_run

#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-185338.4873197' #nodev
#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779' #deven
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-12-193933.6112353'

if len(sys.argv) > 1:
    rundir = sys.argv[1]

zero_ground = False
measure_deviation = False
use_volts = False
do_cut = False
black_dots = False

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

plt.rcParams['font.size'] = 14
title = None
outname = None
make_debug_plots = False
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
    if lrg == ('-s'):
        scale = sys.argv[idx + 1]
        if not scale in scale_map:
            print(f"Error: unrecognized scale {scale}")
            sys.exit(111)
        idx += 1
    if lrg == '-b': black_dots = True
    if lrg == '-t':
        title = sys.argv[idx + 1]
        idx += 1
    if lrg == '-o':
        outname = sys.argv[idx + 1]
        idx += 1
    if lrg == '-d':
        make_debug_plots = True

kV_from_V = 1E-03 # kV/V
scale_mult = scale_map[scale][0]
scale_lab = scale_map[scale][1]


run = aef_run.aef_run(rundir)
run_str = run.run

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]
print(Ez)
print(states)
Ezs = df[Ez]
print(Ezs)
mid_idx = len(Ezs) // 2
Ez_mid = Ezs[mid_idx] * kV_from_V
Ezs_kV = np.array(df[Ez]) * kV_from_V
going_backwards = False

if make_debug_plots:
    os.makedirs(os.path.join(rundir, "debug_plots"), exist_ok=True)

## Dot-product state tracking starts at the E-field
## Furthest away from zero and goes "backwards" towards zero
if abs(Ezs[0]) > abs(Ezs[1]):
    going_backwards = True
    # flip 
    #Ezs = Ezs[::-1]
    #dE_gnds

trans_table = aef_run.state_translation_table(run)

## This now takes the state translation table into account
def get_energies(sdx):
    n_Ezs = len(Ezs)
    Es = np.zeros(n_Ezs)
    for i in range(n_Ezs):
        E_z = Ezs[i]
        edx = trans_table.edx_from_sdx_Ez(E_z, sdx)
        Es[i] = df[f'E{edx}'][i]
    return Es
## 
## suspect MDA is actually backwards
# +Z oriented states start at 0
## expected structure -- doublet + doublet + quadruplet
grp_size = 8
bidx_pz = 0
didx_f1_0 = 0 # delta_idx to reach f_1 = 0
didx_f1_1 = 2 # delta_idx to reach f_1 = 1
didx_f1_1_f_half = 2 # delta idx to reach f_1 = 1, f = 1/2
didx_f1_1_f_thlf = 4 # delta idx to reach f_1 = 1, f = 3/2 -- thlf == 3/2


## Temporarily hardcode things
## Translation table between old and new indicies
## Old index | What it was | new index | what it is now                  |
##     0     | f=0, m_f = 0|     0     |  f_1 = 0, f = 1/2, m_f = -1/2   |
##     1     | f=1, m_f = 0|     2     |  f_1 = 1, f = 1/2, m_f = -1/2   |
##     2     | f=1, m_f =-1|     4     |  f_1 = 1, f = 1/2, m_f = -3/2   |
##     3     | f=1, m_f = 1|     7     |  f_1 = 1, f = 1/2, m_f = +3/2   |   

Ep0s = get_energies(bidx_pz + 0) # doublet f_1 = 0, f = 1/2
Ep1s = get_energies(bidx_pz + 1)
Ep2s = get_energies(bidx_pz + 2) # doublet f_1 = 1, f = 1/2
Ep3s = get_energies(bidx_pz + 3)
Ep4s = get_energies(bidx_pz + 4) # quadruplet f_1 = 1, f = 3/2
Ep5s = get_energies(bidx_pz + 5)
Ep6s = get_energies(bidx_pz + 6)
Ep7s = get_energies(bidx_pz + 7)

if make_debug_plots:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title("Debug plot: Whole Bottom group energies")
    Ess = [get_energies(bidx_pz + idx) for idx in range(6*grp_size)]
    idx = 0
    for Es in Ess:
        plt.plot(Ezs_kV, Es * scale_mult, label=f'State Index {idx}')
        idx += 1
    plt.legend()
    plt.ylabel(f"Energy ({scale_lab})")
    plt.xlabel(f"Externally-Applied Electric Field (kV/cm)")
    plt.savefig(os.path.join(rundir, "debug_plots", "debug_plot_0_bottom_group_spect.png"))
    plt.savefig(os.path.join(rundir, "debug_plots", "debug_plot_0_bottom_group_spect.pdf"))
    plt.savefig(os.path.join(rundir, "debug_plots", "debug_plot_0_bottom_group_spect.svg"))
    plt.show()
    sys.exit(0)

Eps = [get_energies(bidx_pz + idx) for idx in range(grp_size)]

Em0s = get_energies(bidx_pz + grp_size)
Em1s = get_energies(bidx_pz + grp_size + 1)
# -z oriented states
bidx_nz = grp_size * 5 #20 # 4 (for -z) + 4*4 (ffor +-x, +-y)
Ens = [get_energies(bidx_nz + idx) for idx in range(grp_size)]
En0s = get_energies(bidx_nz + 0)
En1s = get_energies(bidx_nz + 1)
En2s = get_energies(bidx_nz + 2)
En3s = get_energies(bidx_nz + 3)
En4s = get_energies(bidx_nz + 4)
En5s = get_energies(bidx_nz + 5)
En6s = get_energies(bidx_nz + 6)
En7s = get_energies(bidx_nz + 7)

delta_10_pmf1 = Ep4s[1] - Ep0s[1] # was 2,0
delta_10_pmft = Ep7s[1] - Ep0s[1] # was 3,0
delta_10_nmf1 = En4s[1] - En0s[1] # was 2,0
delta_10_nmft = En7s[1] - En0s[1] # was 3,0


dn0 = Ep2s[1] - Ep0s[1]#np.average(En1s - En0s)
dp0 = En2s[1] - En0s[1]#np.average(Ep1s - Ep0s)

### Make equivalent to PRA fig 3
props = dict(boxstyle='round,pad=0.2', facecolor='wheat', alpha=0.5)
x_pix = 698
y_pix = 1025
fig,ax = plt.subplots(4, 1, figsize=(x_pix / 100.0, y_pix / 100.0), sharex = True, sharey = False)
# part a
gca = plt.subplot(4, 1, 1)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'Stark Shift of Lowest-energy group of states'
color = []
for i in range(4): color.append('b')
for i in range(16): color.append('g')
for i in range(4): color.append('r')
gca.text(0.015, 0.12, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)

#df.plot(Ez[:24], states[:24], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = False)
plt.ylabel(f'Energy ({scale_lab})')
#plt.plot(df[Ez][1:], df[states[:24]][1:] * scale_mult)
plt.plot(Ezs_kV, Ep0s * scale_mult, color='b')
plt.annotate('$+\hat{Z}$', xy=(Ez_mid, Ep0s[mid_idx - 1] * scale_mult), xycoords='data', xytext=(1.5, 1.5), color='b', textcoords='offset points')
plt.plot(Ezs_kV, Em0s * scale_mult, color='g')
plt.annotate('$+\hat{X},-\hat{X},+\hat{Y},-\hat{Y}$', xy=(Ez_mid, Em0s[mid_idx - 1] * scale_mult), xycoords='data', xytext=(1.5, 2.5), color='g', textcoords='offset points')
plt.plot(Ezs_kV, En0s * scale_mult, color='r')
plt.annotate('$-\hat{Z}$', xy=(Ez_mid, En0s[mid_idx - 1] * scale_mult), xycoords='data', xytext=(1.5, 5.5), color='r', textcoords='offset points')

#for i,j in enumerate(gca.lines):
#    j.set_color(color[i])

# part b
gca = plt.subplot(4, 1, 2)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'Hyperfine shift of $f_1=1$ above $f_1=0$\nfor +Z part of the lowest-energy group'
gca.text(0.05, 0.25, textstr, transform=gca.transAxes, fontsize=12, verticalalignment='top', bbox=props)
plt.ylabel('[Energy ($f_1$=1) - Energy ($f_1$=0)](MHz)')
s = '$f_1 = 1, f = 1/2$'
plt.plot(Ezs_kV[1:], (En2s-En0s)[1:], label=s)
plt.annotate(s, xy=(Ez_mid, (En2s-En0s)[mid_idx - 1]), xycoords='data', xytext=(1.5, 5.5), color='k', textcoords='offset points')
s = '$f_1 = 1, f = 3/2$'
plt.plot(Ezs_kV[1:], (En4s-En0s)[1:], label=s)
plt.annotate(s, xy=(Ez_mid, (En4s-En0s)[mid_idx - 1]), xycoords='data', xytext=(1.5, -9.5), color='k', textcoords='offset points')
plt.plot(Ezs_kV[1:], (En5s-En0s)[1:], label='$f_1=1,f=3/2,m_f=-1/2$')
#plt.legend()
gca.set_ylim([60, 68])

# part c
gca = plt.subplot(4, 1, 3)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'$f_1=1,f=3.2$'
plt.ylabel('[E($f_1$=1)-E($f_1$=0)-$\Delta_{10}(0)$] (kHz)')
gca.text(0.05, 0.15, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
pz_avg = ((Ep4s-Ep0s-delta_10_pmf1) + (Ep7s-Ep0s-delta_10_pmft)) / 2.0
nz_avg = ((En4s-En0s-delta_10_nmf1) + (En7s-En0s-delta_10_nmft)) / 2.0
plt.plot(Ezs_kV[1:], pz_avg[1:]*1000, 'b-', label='+Z')
plt.annotate('+Z', xy=(Ez_mid, pz_avg[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 5.5), color='b', textcoords='offset points')
plt.plot(Ezs_kV[1:], nz_avg[1:]*1000, 'r-', label='-Z')
plt.annotate('-Z', xy=(Ez_mid, nz_avg[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 1.5), color='r', textcoords='offset points')
plt.legend()

# part d
gca = plt.subplot(4, 1, 4)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'$f=1,m_f=0$'
plt.ylabel('[E($f_1$=1)-E($f_1$=0)-$\Delta_{10}(0)$] (kHz)')
gca.text(0.05, 0.15, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
pz = (Ep2s-Ep0s-dn0)  # was Ep1-En0
nz = (En2s-En0s-dp0)  # was En1-En0
plt.plot(Ezs_kV[1:], pz[1:]*1000, 'b-', label='+Z')
plt.annotate('+Z', xy=(Ez_mid, pz[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 1.5), color='b', textcoords='offset points')
plt.plot(Ezs_kV[1:], nz[1:]*1000, 'r-', label='-Z')
plt.annotate('-Z', xy=(Ez_mid, nz[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 5.5), color='r', textcoords='offset points')
plt.legend()

fig.suptitle(f"N=0, $F_1$=0,1 Stark shift for run {run_str}", y=0.999)
plt.subplots_adjust(bottom=0.05, right=0.990, top=0.97, left = 0.114, hspace = 0.0)
plt.xlabel("Externally-Applied Electric field Strength (kV/cm)")

plt.savefig(os.path.join(rundir, 'doe_plot_pra_aspect_deven.png'))
plt.show()