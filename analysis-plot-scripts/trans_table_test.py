#!/usr/bin/env python3
## trans_table_test.py -- intended for testing the state_translation_table
## features of the aef_run python module
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

trans_table = aef_run.state_translation_table(run)
print("Finished parsing state translation table")
## expected structure -- doublet (f_1 = 1, f = 1/2) + quadruplet (f_1 = 1, f = 3/2) + doublet (f_1 = 0, f = 1/2)
## because B_F is negative
grp_size = 8
bidx_pz = 0
didx_f1_1 = 0 # delta_idx to reach f_1 = 1
didx_f1_1_f_half = 0 # delta idx to reach f_1 = 1, f = 1/2
didx_f1_1_f_thlf = 2 # delta idx to reach f_1 = 1, f = 3/2 -- thlf == 3/2
didx_f1_0 = 6 # delta_idx to reach f_1 = 0

plt.figure(figsize=(13.66,9.00))
for sdx in range(10):
    plt.scatter(trans_table.E_z_list / 1000, trans_table.edx_from_sdx_arr[:,sdx], label=f"Invariant State Index {sdx}")
plt.title("Translation Table")
plt.ylabel("Energy Eigenstate Index")
plt.xlabel("Externally-Applied Electric Field (kV/cm)")
plt.show()