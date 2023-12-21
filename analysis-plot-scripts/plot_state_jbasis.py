#!/usr/bin/env python3
## plot_state_jbasis.py -- plots the probability of each j-basis ket in
## a given state. 
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
import math
import cmath
from random import triangular
import sys
import os
import os.path
import re
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import numpy.linalg as npla
import pandas as pd
import numba

import aef_run
import baf_state
import triangular_state_plotter


#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-15-195553.6276441'
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-08-09-135429.5590762'
if len(sys.argv) > 1:
    rundir = sys.argv[1]
run = aef_run.aef_run(rundir)

#
coeffdir = run.get_coeff_dir()

Ez = 500 
#Ez = 50000
df = run.parse_state_coeffs(Ez)#50000)
#df = pd.DataFrame()
pltdir = os.path.join(coeffdir, f'{Ez}')
os.makedirs(pltdir, exist_ok=True)

nBasisElts = (len(df.keys()) - 1) // 2

hsts = baf_state.make_hyperfine_states(run.nmax)

@numba.njit
def mag_sq(vec):
    vabs = np.absolute(vec)
    return vabs*vabs

def invert_qsq(expect_qsq):
    return (np.sqrt(4 * expect_qsq + 1.0) - 1.0) / 2.0;

def expect_qsq(vec):
    hst = baf_state.HyperfineState(n=0,j=0,f=0,m_f=0)
    prob_vec = mag_sq(vec)
    prob_tot = 0
    for kidx in range(len(vec)):
        ket:baf_state.HyperfineState = hsts[kidx]
        prob_tot += prob_vec[kidx] ## delibrately not using sum
        hst.n += prob_vec[kidx] * ket.n * (ket.n + 1)
        hst.j += prob_vec[kidx] * ket.j * (ket.j + 1)
        hst.f += prob_vec[kidx] * ket.f * (ket.f + 1)
        hst.m_f += prob_vec[kidx] * ket.m_f
    hst.n = invert_qsq(hst.n / prob_tot)
    hst.j = invert_qsq(hst.j / prob_tot)
    hst.f = invert_qsq(hst.f / prob_tot)
    hst.m_f /= prob_tot
    return hst

def get_names(kidx):
    return (f'Re(<j_{kidx}|E_n>)', f'Im(<j_{kidx}|E_n>')

def read_state(df, idx):
    ket = np.zeros(nBasisElts, dtype=np.complex128)
    row = (df.T[idx])
    print(row)
    for kidx in range(nBasisElts):
        names = get_names(kidx)
        ket[kidx] = row[names[0]] + 1j*row[names[1]]
    prob_tot = sum(mag_sq(ket))
    dev = 1 - prob_tot
    rel_dev = dev / prob_tot
    print(f'magsq of {idx} is {prob_tot} w/ dev {dev}, rel_dev {rel_dev}')
    return ket

def plot_state(st, stnam, fig:plt.Figure=None, typ = 'mag', cmap='viridis'):
    if fig == None:
        fig = plt.figure(figsize=(19.2, 12.42))#11.0*19.2/17.0))#19.2, 16.8)) #fig = plt.gcf()
    # turn argument into colormap
    norm = colors.LogNorm(1e-19, 1)
    cmap = mpl.colormaps.get_cmap(cmap)
    def f_mag(ax:plt.Axes, hst:baf_state.HyperfineState, idx:int, x:float, y:float):
        k = npla.norm(st[idx])**2
        #if k == 0: k = 1e-38
        # get color from map
        color = cmap(norm(k))
        print(idx, hst, k, norm(k))
        ax.annotate(f'{k:.2e}\n{idx}', (x, y - 0.3), ha='center', size=8)
        ax.plot(x, y, 'o', color=color, markersize=12)
        return 1
    njs = 0
    for n in range(4):
        nl = triangular_state_plotter.nlevel(n)
        nl.draw_boxes(plt.gca(), njs, f_mag)
        prev_njs = njs
        njs += 5 if (n > 0) else 3
        #plt.hlines(-5, prev_njs, njs - 0.5, color='k', linestyle='-', label=f'n={n}')
        dx = ((njs - 0.5) - (prev_njs)) / 2
        x = prev_njs + 0.5
        #plt.text(x + dx, -4.5, f'n={n}', ha='right', va='center')
    plt.title(f'J-basis-ket state plot for state {stnam} of run {run.run}, E_z = {Ez} V/cm')
    plt.ylabel('m_f')
    plt.xlabel('Arbitrary')
    print("state coeff 0: ", st[0], ' sq val: ', npla.norm(st[0])**2)
    print("state coeff 6: ", st[6], ' sq val: ', npla.norm(st[6])**2)
    plt.colorbar(mplcm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.gca())
    onam = os.path.join(run.cdir, f'{Ez}', f'{stnam}.png')
    plt.savefig(onam)#'state_j-basis-plot.png')

if __name__ == '__main__':
    print(f'Plotting states from run {run.run}')

    pz_f00 = read_state(df, 0)
    pz_f10 = read_state(df, 1)
    pz_f1t = read_state(df, 2)
    pz_f11 = read_state(df, 3)

    nz_f00 = read_state(df, 4)#20)
    nz_f10 = read_state(df, 5)#21)
    nz_f1t = read_state(df, 6)#23)
    nz_f11 = read_state(df, 7)#22)

    print('Read states, starting plotting')
    fig = None
    #fig = plt.figure(figsize=(19.2, 16.8))
    plot_state(pz_f00, 'pz_f00', fig)
    plot_state(pz_f10, 'pz_f10', fig)
    plot_state(pz_f1t, 'pz_f1t', fig)
    plot_state(pz_f11, 'pz_f11', fig)
    #nz
    plot_state(nz_f00, 'nz_f00', fig)
    plot_state(nz_f10, 'nz_f10', fig)
    plot_state(nz_f1t, 'nz_f1t', fig)
    plot_state(nz_f11, 'nz_f11', fig)
    plt.show()
    