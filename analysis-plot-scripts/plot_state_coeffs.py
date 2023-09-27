
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

rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-15-195553.6276441'
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

def get_names(kidx):
    return (f'Re(<j_{kidx}|E_n>)', f'Im(<j_{kidx}|E_n>')

def read_state(df, idx):
    ket = np.zeros(nBasisElts, dtype=np.complex128)
    row = (df.T[idx])
    print(row)
    for kidx in range(nBasisElts):
        names = get_names(kidx)
        ket[kidx] = row[names[0]] + 1j*row[names[1]]
    return ket

pz_f00 = read_state(df, 0)
pz_f10 = read_state(df, 1)
pz_f1t = read_state(df, 2)
pz_f11 = read_state(df, 3)

nz_f00 = read_state(df, 4)#20)
nz_f10 = read_state(df, 5)#21)
nz_f1t = read_state(df, 6)#23)
nz_f11 = read_state(df, 7)#22)


if 0:
    ## real part
    fig = plt.figure(figsize=(12.00, 9.00))
    plt.plot(np.real(pz_f00), label='F00')
    plt.plot(np.real(pz_f10), label='F10')
    plt.plot(np.real(pz_f11), label='F11')
    plt.plot(np.real(pz_f1t), label='F1T')
    plt.legend()
    plt.show()


## magnitude
fig = plt.figure(figsize=(12.00, 9.00))
plt.plot(np.abs(pz_f00)**2, label='0,+Z,F00')
plt.plot(np.abs(pz_f10)**2, label='0,+Z,F10')
plt.plot(np.abs(pz_f11)**2, label='0,+Z,F11')
plt.plot(np.abs(pz_f1t)**2, label='0,+Z,F1T')
plt.yscale('log')
plt.title('+Z oriented Coefficent Magnitude-Squared vs Basis Index')
plt.xlabel('J-Basis Ket index')
plt.ylabel(r'$|\left<j_i|E_n\right>|^2$')
plt.legend()
plt.savefig(os.path.join(pltdir, 'bottom_grp_pz_mag.png'))
#plt.show()

fig = plt.figure(figsize=(12.00, 9.00))
plt.plot(np.abs(nz_f00)**2, label='0,-Z,F00')
plt.plot(np.abs(nz_f10)**2, label='0,-Z,F10')
plt.plot(np.abs(nz_f11)**2, label='0,-Z,F11')
plt.plot(np.abs(nz_f1t)**2, label='0,-Z,F1T')
plt.yscale('log')
plt.title('-Z oriented Coefficent Magnitude-Squared vs Basis Index')
plt.xlabel('J-Basis Ket index')
plt.ylabel(r'$|\left<j_i|E_n\right>|^2$')
plt.legend()
plt.savefig(os.path.join(pltdir, 'bottom_grp_nz_mag.png'))

def assemble_matrix(bidx):
    mat = np.zeros((nBasisElts, 4), dtype=np.complex128)
    mat[:, 0] = read_state(df, bidx+0) # f00
    mat[:, 1] = read_state(df, bidx+2) # f1t
    mat[:, 2] = read_state(df, bidx+1) # f10
    mat[:, 3] = read_state(df, bidx+3) # f11
    return mat

def prune_matrix(mat, limit=1e-18):
    mdx = 0
    for j in range(mat.shape[1]):
        for i in range(mat.shape[0]):
            if (np.abs(mat[i, j])**2) > limit:
                if i > mdx: mdx = i
            else: mat[i, j] = 0
    return mdx

def plot_line(mat, nam, tnam=None, mdx=None):
    if tnam == None: tnam = nam
    if mdx == None: mdx = mat.shape[0] - 1
    fig = plt.figure(figsize=(12.00, 9.00))
    plt.plot(np.abs(mat[0:mdx,0])**2, 'b-', label=f'{nam},F00')
    plt.plot(np.abs(mat[0:mdx,1])**2, '-', color='orange', label=f'{nam},F1T')
    plt.plot(np.abs(mat[0:mdx,2])**2, 'g-', label=f'{nam},F10')
    plt.plot(np.abs(mat[0:mdx,3])**2, 'r-', label=f'{nam},F11')
    plt.yscale('log')
    plt.title(f'{tnam} \"Pruned\" Coefficent Magnitude-Squared vs Basis Index')
    plt.xlabel('J-Basis Ket index')
    plt.ylabel(r'$|\left<j_i|E_n\right>|^2$')
    plt.legend()
    plt.savefig(os.path.join(pltdir, f'bottom_grp_{nam}_mag.png'))

    
fig = plt.figure(figsize=(12.00, 9.00))
pmat = np.zeros((nBasisElts, 4), dtype=np.complex128)
pmat[:, 0] = pz_f00
pmat[:, 1] = pz_f10
pmat[:, 2] = pz_f1t
pmat[:, 3] = pz_f11
mdx = prune_matrix(pmat, 1e-18)
print(mdx)

plot_line(pmat, '0,+Z', '+Z oriented', mdx)

plt.matshow(np.abs(pmat[0:mdx,:].T)**2, fig, norm='log', aspect='auto')
plt.title('\"Pruned\" Matrix representation of Basis Amplitudes-Squared for Bottom-group +Z states')
plt.xlabel('J-basis Index')
plt.ylabel('Energy eigenstate number')
plt.colorbar()
plt.savefig(os.path.join(pltdir, 'matshow_bg_pz.png'))

fig = plt.figure(figsize=(12.00, 9.00))
nmat = np.zeros((nBasisElts, 4), dtype=np.complex128)
nmat[:, 0] = nz_f00
nmat[:, 1] = nz_f10
nmat[:, 2] = nz_f1t
nmat[:, 3] = nz_f11
mdx = prune_matrix(nmat, 1e-18)
plot_line(nmat, '0,-Z', '-Z oriented', mdx)
plt.matshow(np.abs(nmat[0:mdx,:].T)**2, fig, norm='log', aspect='auto')
plt.title('\"Pruned\" Matrix representation of Basis Amplitudes-Squared for Bottom-group -Z states')
plt.xlabel('J-basis Index')
plt.ylabel('Energy eigenstate number')
plt.colorbar()
plt.savefig(os.path.join(pltdir, 'matshow_bg_nz.png'))

plt.show()