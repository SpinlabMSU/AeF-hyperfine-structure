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


rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-18-190613.6465128'
run = os.path.split(rundir)[1]
if len(sys.argv) > 1:
    rundir = sys.argv[1]

starkpath = os.path.join(rundir, 'stark_shift_gnd.csv')
df = pd.read_csv(starkpath)

## csv keys -- note that the spaces in front are intentional
key_E = 'E-field (V/cm)'
key_dEgnd = ' dE_gnd'
key_dEf1t = ' dE_f1t'
key_dEf10 = ' dE_f10'
key_dEf11 = ' dE_f11'

fig = plt.figure(figsize=(13.66, 9.00))
gca = plt.subplot(2, 1, 1)
gca.set_title(f'F=1 Ground-relative Stark shift for run {run}')
df.plot(key_E, [key_dEf1t, key_dEf10, key_dEf11], ax=gca)
gca = plt.subplot(2, 1, 2)
gca.set_title(f'Ground State Stark Shift run {run}')
df.plot(key_E, key_dEgnd, ax=gca, ylabel='Energy (MHz)')
plt.show()