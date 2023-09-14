#!/usr/bin/env python3
### collate the runtimes of a series of runs

# system imports
import math
import cmath
import sys
import os
import os.path
import re
# scipy
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npla
import pandas as pd
import numba

# user-defined modules
import aef_run

indir = r'C:\Users\nusgart\source\repos\SpinlabMSU\AeF-hyperfine-structure\output'
if len(sys.argv) > 1:
    indir = sys.argv[1]

runlist = aef_run.find_runs(indir)

# matrix element calc times
mtimes_deven = {}
mtimes_nodev = {}

# stark loop times
stimes_deven = {}
stimes_nodev = {}

for run in runlist:
    if not run.valid:
        print(f"run {run.run} is not valid, please check the following path: {run.path}", file=sys.stderr)
    nmax = run.nmax
    if run.dev_en:
        mtimes_deven[nmax] = run.mat_elt_dur
        stimes_deven[nmax] = run.stk_lop_dur
    else:
        mtimes_nodev[nmax] = run.mat_elt_dur
        stimes_nodev[nmax] = run.stk_lop_dur

nmaxes = set(mtimes_deven.keys()).union(mtimes_nodev.keys())

## config var
tabularx = True

print(r'\begin{table}[]')
print('\t\\centering')
if not tabularx:
    print('\t\\begin{tabular}{|c|c|c|c|c|}\\hline')
else:
    print(r'''\begin{tabularx}{1.0\textwidth}{ 
  | >{\centering\arraybackslash}X 
  | >{\centering\arraybackslash}X
  | >{\centering\arraybackslash}X
  | >{\centering\arraybackslash}X
  | >{\centering\arraybackslash}X|}\hline%{|c|c|c|c|c|}\hline''')
print('\t\tMaximum n value & Vacuum matrix element calculation duration (s)& Vacuum Stark loop duration (s) '+
      '& In-matrix matrix element calculation duration (s) & In-matrix stark loop duration (s)\\\\\\hline')
for nmax in nmaxes:
    mtime_nodev = -1
    stime_nodev = -1
    mtime_deven = -1
    stime_deven = -1

    if nmax in mtimes_nodev:
        mtime_nodev = mtimes_nodev[nmax]
        stime_nodev = stimes_nodev[nmax]
    #else: print('no nmax for nodev')
    if nmax in mtimes_deven:
        mtime_deven = mtimes_deven[nmax]
        stime_deven = stimes_deven[nmax]
    #else: print('no nmax for deven')
    #
    dtab = '\t\t'
    print(fr'{dtab}{nmax}&{mtime_nodev}&{stime_nodev}&{mtime_deven}&{stime_deven}\\\hline')
print('''
    \end{tabular}
    \caption{Sample runtimes for a variety of values of $n_{max}$ with Devonshire enabled and disabled.}
    \label{tab:runtimes}
\end{table}''')
#
############ CUT HERE #############
sys.exit(0)
## old version
#!/usr/bin/env python3
import os
import sys
import os.path
from dataclasses import dataclass
import re

def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

if True:
  marker = 'have taken '
  mlen = len(marker)
  def get_duration_from_line(lline):
    pdx = lline.index(marker) + mlen
    pstr = lline[pdx:]
    dur_str = pstr.split(' ')[0]
    return float(dur_str)

def get_ssv_val(line, marker, typ=int):
  "parse a space-seperated value line, returning the value preceeded by marker"
  pdx = line.index(marker) + len(marker) + 1
  line_v_to_eol = line[pdx:]
  val_str = re.split(r'[,\s]', line_v_to_eol)[0]
  return typ(val_str)

@dataclass
class RunInfo:
  run:str
  nmax: int
  mat_elt_dur: float
  stk_lop_dur: float
  
  def __getitem__(self, key):
    return getattr(self, key)
  def __setitem__(self, key, value):
    setattr(self, key, value)

def process_dir(odir):
  run = os.path.basename(odir)
  olog = open(os.path.join(odir, 'out.log'))
  
  ifo = RunInfo(run, -1, -1, -1)
  
  for line in olog:
    lline = line.lower()
    if lline.startswith('finished matrix elt'):
      # mat elt completion line --> parse for duration
      ifo.mat_elt_dur = get_duration_from_line(lline)
    elif lline.startswith('completed stark loop'):
      # stark loop completion line --> parse for duration
      ifo.stk_lop_dur = get_duration_from_line(lline)
    elif lline.startswith('nmax is '):
      ifo.nmax = get_ssv_val(lline, 'nmax is', int)
    else: pass
   
  return ifo
 
def process_dirs(base_dir):
  ifos = {}
  for entry in os.scandir(base_dir):
    if not entry.is_dir(): continue
    print(f"[dbg] processing {entry.name}")
    ifo = process_dir(entry.path)
    ifos[ifo.run] = ifo
    #print(f"[dbg] processed {entry.name} as ")
  return ifos
print (process_dir('./2023-08-27-113340.1263312'))

data = process_dirs('.')

print(data)

print('Run,nmax,Matrix Element Calculation duration (s),Stark Loop duration (s)')
for run in data.keys():
  rundat = data[run]
  print(f"{run},{rundat.nmax},{rundat.mat_elt_dur},{rundat.stk_lop_dur}")

print()
print()
print()
#f = open('', 'w')
#print('
for run in data.keys():
  rundat = data[run]
  latex_eol = r'\\\hline'
  print(f"{rundat.nmax}&{rundat.mat_elt_dur}&{rundat.stk_lop_dur}{latex_eol}")