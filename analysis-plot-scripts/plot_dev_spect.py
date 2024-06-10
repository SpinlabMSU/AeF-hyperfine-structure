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


#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-185338.4873197' #nodev
#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779' #deven
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-12-193933.6112353'

if len(sys.argv) > 1:
    rundir = sys.argv[1]
run = os.path.split(rundir)[1]

starkpath = os.path.join(rundir, 'stark_spectrum.csv')
df = pd.read_csv(starkpath)

Ez = df.keys()[1]
states = df.keys()[2:]
print(Ez)
print(states)
Ezs = df[Ez]
mid_idx = len(Ezs) // 2
Ez_mid = Ezs[mid_idx]

do_extras = False
# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Energy Spectrum for run {run}")
    df.plot(Ez, states, ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
    #plt.show()
    plt.close(fig)

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Energy Spectrum for run {run}")
    df.plot(Ez[:24], states[:24], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = True)
    plt.savefig(os.path.join(rundir, 'spectrum_bottom_group.png'))
    #plt.show()
    plt.close(fig)

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

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Descending part for run {run}")
    plt.plot(df[Ez][1:], (Ep1s-Ep0s)[1:], label='First Excited state')
    plt.plot(df[Ez][1:], (Ep2s-Ep0s)[1:], label='Second Excited state')
    plt.plot(df[Ez][1:], (Ep3s-Ep0s)[1:], label='Third Excited state')
    plt.legend()
    plt.ylabel('Energy (MHz)')
    plt.xlabel("Electric Field (V/cm)")
    #df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_bottom_negz.png'))
    plt.close(fig)
    #plt.show()

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Energy of bottom ascending part for run {run}")
    plt.plot(df[Ez][1:], (En1s-En0s)[1:], label='First Excited state')
    plt.plot(df[Ez][1:], (En2s-En0s)[1:], label='Second Excited state')
    plt.plot(df[Ez][1:], (En3s-En0s)[1:], label='Third Excited state')
    plt.legend()
    plt.ylabel('Energy (MHz)')
    plt.xlabel("Electric Field (V/cm)")
    #df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_bottom_posz.png'))
    plt.close(fig)


delta_10_mf0 = Ep1s[1] - Ep0s[1]
delta_10_mf1 = Ep2s[1] - Ep0s[1]
delta_10_mft = Ep3s[1] - Ep0s[1]

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Descending part for run {run}")
    plt.plot(df[Ez][1:], (Ep1s-Ep0s-delta_10_mf0)[1:]*1000, label='First Excited state')
    plt.plot(df[Ez][1:], (Ep2s-Ep0s-delta_10_mf1)[1:]*1000, label='Second Excited state')
    plt.plot(df[Ez][1:], (Ep3s-Ep0s-delta_10_mft)[1:]*1000, label='Third Excited state')
    plt.legend()
    plt.ylabel('Energy (kHz)')
    plt.xlabel("Electric Field (V/cm)")
    #df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_bottom_stks_d10.png'))
    plt.close(fig)
    #plt.show()

delta_10_pmf1 = Ep2s[1] - Ep0s[1]
delta_10_pmft = Ep3s[1] - Ep0s[1]
delta_10_nmf1 = En2s[1] - En0s[1]
delta_10_nmft = En3s[1] - En0s[1]

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Equiv of EDM3 RMP Fig 3c for run {run}")
    plt.plot(df[Ez][1:], (Ep2s-Ep0s-delta_10_pmf1)[1:]*1000, label='+Z,f=1,m_f=?1?')
    plt.plot(df[Ez][1:], (Ep3s-Ep0s-delta_10_pmft)[1:]*1000, label='+z,f=1,m_f=?-1?')
    plt.plot(df[Ez][1:], (En2s-En0s-delta_10_nmf1)[1:]*1000, label='-Z,f=1,m_f=?1?')
    plt.plot(df[Ez][1:], (En3s-En0s-delta_10_nmft)[1:]*1000, label='-Z,f=1,m_f=?-1?')
    plt.legend()
    plt.ylabel('Energy (kHz)')
    plt.xlabel("Electric Field (V/cm)")
    #df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_rmp_3c.png'))
    plt.close(fig)


dn0 = Ep1s[1] - Ep0s[1]#np.average(En1s - En0s)
dp0 = En1s[1] - En0s[1]#np.average(Ep1s - Ep0s)

if do_extras:
    fig = plt.figure(figsize=(13.66, 9.00))
    plt.title(f"Equiv to EDM3 RMP Fig3d for run {run}")
    plt.plot(df[Ez][1:], (Ep1s-Ep0s-dn0)[1:]*1000, label='+Z,f=1,m_f=0')
    plt.plot(df[Ez][1:], (En1s-En0s-dp0)[1:]*1000, label='-Z,f=1,m_f=0')
    plt.legend()
    plt.ylabel('Energy (kHz)')
    plt.xlabel("Electric Field (V/cm)")
    #df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
    plt.savefig(os.path.join(rundir, 'spectrum_bottom_fig3d.png'))
    #plt.show()
    plt.close(fig)

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
plt.ylabel('Energy (GHz)')
#plt.plot(df[Ez][1:], df[states[:24]][1:] / 1000)
plt.plot(df[Ez], Ep0s/1000, color='b')
plt.annotate('$+\hat{Z}$', xy=(Ez_mid, Ep0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 1.5), color='b', textcoords='offset points')
plt.plot(df[Ez], Em0s/1000, color='g')
plt.annotate('$+\hat{X},-\hat{X},+\hat{Y},-\hat{Y}$', xy=(Ez_mid, Em0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 2.5), color='g', textcoords='offset points')
plt.plot(df[Ez], En0s/1000, color='r')
plt.annotate('$-\hat{Z}$', xy=(Ez_mid, En0s[mid_idx - 1]/1000), xycoords='data', xytext=(1.5, 5.5), color='r', textcoords='offset points')

#for i,j in enumerate(gca.lines):
#    j.set_color(color[i])

# part b
gca = plt.subplot(4, 1, 2)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'Hyperfine shift of f=1 above f=0 for +Z part of the lowest-energy group'
gca.text(0.01, 0.15, textstr, transform=gca.transAxes, fontsize=12, verticalalignment='top', bbox=props)
plt.ylabel('[Energy (f=1) - Energy (f=0)] (MHz)')
plt.plot(df[Ez][1:], (En1s-En0s)[1:], label='$f=1,m_f=0$')
plt.annotate('$f=1,m_f=0$', xy=(Ez_mid, (En1s-En0s)[mid_idx - 1]), xycoords='data', xytext=(1.5, 5.5), color='k', textcoords='offset points')
plt.plot(df[Ez][1:], (En2s-En0s)[1:], label='$f=1,m_f=1$')
plt.annotate('$f=1,m_f=\\pm1$', xy=(Ez_mid, (En2s-En0s)[mid_idx - 1]), xycoords='data', xytext=(1.5, -9.5), color='k', textcoords='offset points')
plt.plot(df[Ez][1:], (En3s-En0s)[1:], label='$f=1,m_f=-1$')
#plt.legend()
gca.set_ylim([60, 68])

# part c
gca = plt.subplot(4, 1, 3)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'$f=1,m_f=\\pm1$'
plt.ylabel('[E(f=1)-E(f=0)-$\Delta_{10}(0)$] (kHz)')
gca.text(0.05, 0.15, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
pz_avg = ((Ep2s-Ep0s-delta_10_pmf1) + (Ep3s-Ep0s-delta_10_pmft)) / 2.0
nz_avg = ((En2s-En0s-delta_10_nmf1) + (En3s-En0s-delta_10_nmft)) / 2.0
plt.plot(df[Ez][1:], pz_avg[1:]*1000, 'b-', label='+Z')
plt.annotate('+Z', xy=(Ez_mid, pz_avg[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 5.5), color='b', textcoords='offset points')
plt.plot(df[Ez][1:], nz_avg[1:]*1000, 'r-', label='-Z')
plt.annotate('-Z', xy=(Ez_mid, nz_avg[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 1.5), color='r', textcoords='offset points')
plt.legend()

# part d
gca = plt.subplot(4, 1, 4)
gca.tick_params(axis='both', which='both', direction='inout')
textstr = f'$f=1,m_f=0$'
plt.ylabel('[E(f=1)-E(f=0)-$\Delta_{10}(0)$] (kHz)')
gca.text(0.05, 0.15, textstr, transform=gca.transAxes, fontsize=14, verticalalignment='top', bbox=props)
pz = (Ep1s-Ep0s-dn0)
nz = (En1s-En0s-dp0)
plt.plot(df[Ez][1:], pz[1:]*1000, 'b-', label='+Z')
plt.annotate('+Z', xy=(Ez_mid, pz[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 1.5), color='b', textcoords='offset points')
plt.plot(df[Ez][1:], nz[1:]*1000, 'r-', label='-Z')
plt.annotate('-Z', xy=(Ez_mid, nz[mid_idx - 1] * 1000), xycoords='data', xytext=(1.5, 5.5), color='r', textcoords='offset points')
plt.legend()

fig.suptitle(f"N=0, F=0,1 Stark shift for run {run}", y=0.999)
plt.subplots_adjust(bottom=0.05, right=0.990, top=0.97, left = 0.114, hspace = 0.0)
plt.xlabel("Externally-Applied Electric field Strength (V/cm)")

plt.savefig(os.path.join(rundir, 'pra_aspect_deven.png'))
plt.show()