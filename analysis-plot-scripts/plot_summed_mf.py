#!/usr/bin/env python3
## plot_state_jbasis.py -- plots the summed probabilities of all of the basis kets
## of a given n,j,f triplet
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
coeffdir = run.get_coeff_dir()

Ez = 500 
#Ez = 50000

if len(sys.argv) > 2:
    try:
        Ez = int(sys.argv[2])
    except:
        Ez = float(sys.argv[2])

df = run.parse_state_coeffs(Ez)
pltdir = os.path.join(coeffdir, f'{Ez}', 'mf_sum')
os.makedirs(pltdir, exist_ok=True)

nBasisElts = (len(df.keys()) - 1) // 2

hsts = baf_state.make_hyperfine_states(run.nmax)

plot_nmax = 4
if len(sys.argv) > 3:
    plot_nmax = int(sys.argv[3])

plot_nmax = run.nmax
no_text = False
if len(sys.argv) > 4:
   no_text = sys.argv[4].lower().startswith('no_text') 

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

def make_njf_list(n_max:int):
    njfs = []    
    for n in range(n_max):
        for j in (n - 0.5, n + 0.5):
            if j < 0: continue
            for f in (j - 0.5, j + 0.5):
                if f < 0: continue
                njfs.append((n, j, f))
    return njfs

njfs = make_njf_list(plot_nmax)
n_njfs = len(njfs)

def plot_state(st, stnam, fig:plt.Figure=None, typ = 'mag', cmap='viridis', nmax=4):
    if fig == None:
        fig = plt.figure(figsize=(19.2, 12.42))#11.0*19.2/17.0))#19.2, 16.8)) #fig = plt.gcf()
    #
    hst = baf_state.HyperfineState(n=0,j=0,f=0,m_f=0)
    njf_sums = np.zeros(n_njfs)
    labels = []
    for idx in range(n_njfs):
        hst.n, hst.j, hst.f = njfs[idx]
        sum_prob = 0.0
        for hst.m_f in np.arange(-hst.f, hst.f + 1, 1):
            prob = np.absolute(st[hst.index()])**2
            sum_prob += prob
        njf_sums[idx] = sum_prob
        labels.append(f"|n={hst.n},j={hst.j},f={hst.f}>")
    line, = plt.plot(njf_sums, 'o')
    ax = plt.gca()
    plt.yscale('log')
    ax.set_ylim(bottom=1e-12)
    plt.title(f'Summed-$m_f$ j-basis state plot for state {stnam} of run {run.run}, E_z = {Ez} V/cm\n'
              + r'Plotting magnitude-squared of $\sum_{m_f}\left<n,j,f,m_f|E_{idx}\right>$')
    plt.ylabel('Summed probability')
    plt.xlabel('Arbitrary')
    ### set up hover handler, based heavily on https://stackoverflow.com/a/47166787/9537054
    annot = ax.annotate("", xy=(0,0), xytext=(-20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    # close over 
    def motion_notify_handler(event):
        is_vis = annot.get_visible()
        if event.inaxes != ax:
            return # not us
        cont, det = line.contains(event)
        if not cont:
            if is_vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()
            return
        ind = det['ind']
        x, y = line.get_data()
        annot.xy = (x[ind[0]], y[ind[0]])
        text = "{}, {}".format(" ".join(list(map(str,ind))), 
                           " ".join([labels[n] for n in ind]))
        annot.set_text(text); annot.get_bbox_patch().set_alpha(0.4)
        annot.set_visible(True); fig.canvas.draw_idle()
    fig.canvas.mpl_connect('motion_notify_event', motion_notify_handler)
    print("state coeff 0: ", st[0], ' sq val: ', npla.norm(st[0])**2)
    print("state coeff 6: ", st[6], ' sq val: ', npla.norm(st[6])**2)
    #plt.colorbar(mplcm.ScalarMappable(norm=norm, cmap=cmap), ax=plt.gca())
    onam = os.path.join(pltdir, f'{stnam}.png')
    plt.savefig(onam)#'state_j-basis-plot.png')

if __name__ == '__main__':
    print(f'Plotting sum over m_f of probability of states from run {run.run}')

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
    plot_state(pz_f00, 'pz_f00', fig, nmax=plot_nmax)
    plot_state(pz_f10, 'pz_f10', fig, nmax=plot_nmax)
    plot_state(pz_f1t, 'pz_f1t', fig, nmax=plot_nmax)
    plot_state(pz_f11, 'pz_f11', fig, nmax=plot_nmax)
    #nz
    plot_state(nz_f00, 'nz_f00', fig, nmax=plot_nmax)
    plot_state(nz_f10, 'nz_f10', fig, nmax=plot_nmax)
    plot_state(nz_f1t, 'nz_f1t', fig, nmax=plot_nmax)
    plot_state(nz_f11, 'nz_f11', fig, nmax=plot_nmax)
    plt.show()
    