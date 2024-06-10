#!/usr/bin/env python3
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
import scipy.constants

import aef_run

MHz_per_K = scipy.constants.k / scipy.constants.h * 1e-6
plt.rcParams.update({'font.size': 16})
plt.rc('figure', titlesize=20)
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-08-09-135429.5590762'
if len(sys.argv) > 1:
    rundir = os.path.abspath(sys.argv[1])

run = aef_run.aef_run(rundir)

csvs = run.list_ifo_csvs()


arr_Ez = []
arr_Dz = []

use_zero_field = True
if len(sys.argv) > 2:
    use_zero_field = not (sys.argv[2].lower().startswith('n'))
    print(f"use_zero_field?: {use_zero_field}")

# Perform abs value operation?
doabs = True
if len(sys.argv) > 3:
    doabs = not (sys.argv[3].lower().startswith('n'))

# custom plot format
pltfmt = 'ro'
if len(sys.argv) > 4:
    pltfmt = sys.argv[4]
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

if True:
    ## this regex is intended to parse numbers matching the bnf found below
    s_nrgx = '([+-]?[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?)'
    cplx_rgx = fr'\s*\({s_nrgx}\s*\+\s*i\s*\*{s_nrgx}\s*\)\s*'
    crgx = re.compile(cplx_rgx)
    global parse_cxx_complex
    def parse_cxx_complex(s):
        """
        Parses a complex number in the C++ format matching the following BNF (note that spaces are optional)
            <complex> ::= "(" <real> " + i * " <real> ")"
            <real> ::= <digseq> | <digitseq>.<digitseq> | [+-]<real> | <real>[eE][+-]?<digitseq>
            <digitseq> ::= <digit> | <digit><digitseq>
            <digit> ::= [0-9]
        """
        m = crgx.match(s)
        r = float(m.group(1))
        im = float(m.group(3))
        return complex(r, im)


dev_en = run.dev_en
dev_K = run.dev_K
dirname = run.get_state_ifo_dir()
## todo replace
re_name = re.compile('info_Ez_.*\\.csv')
trtbl = str.maketrans('i','j','() *')
for entry in os.listdir(dirname):
    print(f"Testing entry {entry}")
    if not re_name.match(entry): continue
    fnam = os.path.join(dirname, entry)
    if not os.path.isfile(fnam): continue
    print(f'Accepted entry {entry}')
    stmp = entry.split('_')[2]
    Ez = float(stmp.replace('.csv', ''))
    if Ez == 0 and not use_zero_field: continue
    print(f"Ez is {Ez}")
    dat = pd.read_csv(fnam, encoding='windows-1252')
    key = None
    for k in dat.keys():
        kl = k.lower()
        if '|dz|' in kl and 're' in kl:
            key = k
            break
    if key == None:
        raise KeyError(f'CSV {fnam} missing dz column')
    str_dz = (dat)[key][0]
    dz = float(str_dz) #parse_cxx_complex(str_dz)
    dz = dz.real
    if doabs:
        dz = abs(dz)
    arr_Ez.append(Ez)
    arr_Dz.append(dz)

print(arr_Ez)
print(arr_Dz)

midx = np.argmin(arr_Dz)
print(f"Minimum orientation {arr_Dz[midx]} at electric field {arr_Ez[midx]} index {midx}")

Ez = np.array(arr_Ez)
Dz = np.array(arr_Dz)

fig = plt.figure(figsize=(19.2, 10.8))#13.66, 7.68))
#status_txt = f"enabled, K={dev_K}" if dev_en else "disabled"
status_txt = f"In-Matrix, K={dev_K/MHz_per_K}" if dev_en else "In-Vacuum"
title_text = f"Degree of Molecular Orientation along the externally-applied electric field vs Externally applied electric field strength, {status_txt}"#(with Devonshire {status_txt})\nrun {run.run}"
plt.title(title_text)
plt.xlabel('Externally applied electric field (V/cm)')
plt.ylabel('Degree of Molecular Orientation (unitless, 0 to 1)')
plt.plot(Ez, Dz, pltfmt)
if use_zero_field:
    plt.savefig(os.path.join(rundir, 'orientation_vs_Ez_plot.png'), bbox_inches="tight")
else:
    plt.savefig(os.path.join(rundir, 'orientation_vs_Ez_plot_no_zero_field.png'), bbox_inches="tight")
plt.show()