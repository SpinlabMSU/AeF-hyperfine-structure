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

# Including a legend isn't particularly useful past a certain number of states
# since it runs off the edge of the plot and the colors repeat anyways
use_legend = True
if len(df[Ez]) > 15:
    use_legend = False

fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Energy Spectrum for run {run}")
df.plot(Ez, states, ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
plt.savefig(os.path.join(rundir, 'spectrum_plot.png'))
#plt.show()

fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Energy Spectrum for run {run}")
df.plot(Ez[:24], states[:24], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = True)
plt.savefig(os.path.join(rundir, 'spectrum_bottom_group.png'))
#plt.show()

##
def get_energies(n):
    return df[f'E{n}']
## plot 
#
bidx_nz = 0
En0s = get_energies(bidx_nz + 0)
En1s = get_energies(bidx_nz + 1)
En2s = get_energies(bidx_nz + 2)
En3s = get_energies(bidx_nz + 3)
# +z oriented states
bidx_pz = 20 # 4 (for -z) + 4*4 (ffor +-x, +-y)
Ep0s = get_energies(bidx_nz + 0)
Ep1s = get_energies(bidx_nz + 1)
Ep2s = get_energies(bidx_nz + 2)
Ep3s = get_energies(bidx_nz + 3)


fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Descending part for run {run}")
plt.plot(df[Ez][1:], (En1s-En0s)[1:], label='First Excited state')
plt.plot(df[Ez][1:], (En2s-En0s)[1:], label='Second Excited state')
plt.plot(df[Ez][1:], (En3s-En0s)[1:], label='Third Excited state')
plt.legend()
plt.ylabel('Energy (MHz)')
plt.xlabel("Electric Field (V/cm)")
#df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
plt.savefig(os.path.join(rundir, 'spectrum_bottom_negz.png'))
#plt.show()

delta_10_mf0 = En1s[1] - En0s[1]
delta_10_mf1 = En2s[1] - En0s[1]
delta_10_mft = En3s[1] - En0s[1]

fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Descending part for run {run}")
plt.plot(df[Ez][1:], (En1s-En0s-delta_10_mf0)[1:]*1000, label='First Excited state')
plt.plot(df[Ez][1:], (En2s-En0s-delta_10_mf1)[1:]*1000, label='Second Excited state')
plt.plot(df[Ez][1:], (En3s-En0s-delta_10_mft)[1:]*1000, label='Third Excited state')
plt.legend()
plt.ylabel('Energy (kHz)')
plt.xlabel("Electric Field (V/cm)")
#df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
plt.savefig(os.path.join(rundir, 'spectrum_bottom_stks_d10.png'))
#plt.show()


dn0 = np.average(En1s - En0s)
dp0 = np.average(Ep1s - Ep0s)
fig = plt.figure(figsize=(13.66, 9.00))
plt.title(f"Equiv to EDM3 RMP Fig3d for run {run}")
plt.plot(df[Ez][1:], (En1s-En0s-dn0)[1:]*1000, label='-Z')
plt.plot(df[Ez][1:], (Ep1s-Ep0s-dp0)[1:]*1000, label='+Z')
plt.legend()
plt.ylabel('Energy (kHz)')
plt.xlabel("Electric Field (V/cm)")
#df.plot(Ez, states[:4], ylabel = 'Energy (MHz)', ax=plt.gca(), legend = use_legend)
plt.savefig(os.path.join(rundir, 'spectrum_bottom_fig3d.png'))
plt.show()
