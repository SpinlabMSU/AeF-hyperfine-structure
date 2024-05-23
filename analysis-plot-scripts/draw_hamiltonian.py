#!/usr/bin/env python3
## draw_hamiltonian.py -- draws 
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

#rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-09-15-195553.6276441'
rundir = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-08-09-135429.5590762'
if len(sys.argv) > 1:
    rundir = sys.argv[1]
rundir = os.path.abspath(rundir)

run = aef_run.aef_run(rundir)
mat_file = os.path.join(run.path, 'matrix.dat')

## attempt to read