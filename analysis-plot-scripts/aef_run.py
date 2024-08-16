#!/usr/bin/env python3
## aef_run.py: 
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
import datetime
import re
#
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npla
import pandas as pd
import numba

_run_regex = re.compile(rf'\d+-\d\d?-\d\d?-\d+.\d+')
#numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
def matches_run_spec(dirname):
    return _run_regex.match(dirname)

def get_ssv_val(line, marker, typ=int):
  "parse a space-seperated value line, returning the value preceeded by marker"
  pdx = line.index(marker) + len(marker) + 1
  line_v_to_eol = line[pdx:]
  val_str = re.split(r'[,\s]', line_v_to_eol)[0]
  return typ(val_str)

class aef_run(object):
    
    def __init__(self, dir_path):
        dir_path = dir_path.rstrip('/\\')
        self.path = dir_path
        self.run = os.path.split(dir_path)[1]
        ts = self.run.split('.')
        print(self.run)
        us = ts[1][0:5]
        rs = ts[0] + '.' + us
        self.timestamp = datetime.datetime.strptime(rs, "%Y-%m-%d-%H%M%S.%f")

        self.log_path = os.path.join(self.path, 'out.log')

        self.valid = self.check_valid()
        self.dev_en = False
        self.dev_K = math.nan
        self.nmax = -1 # n is the orbitorotational quantum number
        self.max_E_z = -math.inf
        self.calc_E_z = math.nan

        self.stk_lop_dur = -1
        self.mat_elt_dur = -1

        self.cdir = os.path.join(self.path, 'state_coeffs')
        self.state_info_dir = os.path.join(self.path, 'state_info')

        if not os.path.exists(self.state_info_dir):
            ## old name was devonshire_info because this information was originally only computed when the
            ## devonshire potential was enabled to check for orientation-locking
            self.state_info_dir = os.path.join(self.path, 'devonshire_info')

        if self.valid:
            self.parse_log()

    def check_valid(self):
        run_exists = os.path.exists(self.path) and os.path.isdir(self.path)
        log_exists = os.path.exists(self.log_path) and os.path.isfile(self.log_path)
        self.valid = run_exists and log_exists
        return self.valid

    def parse_log(self):
        f = open(self.log_path, 'r')
        found_param_line = False
        for line in f:
            n_search = "nmax is"
            e_search = "E_z is" # externally applied
            k_search = "K is" # devonshire coupling constant named K, look for enabled
            Emax_search = "Electric field strength is"
            if n_search in line and k_search in line:
                print(f"Found parameter line \"{line}\"")

                #ldx = line.find(n_search) + len(n_search)
                self.nmax = get_ssv_val(line, n_search, int)#self.nmax = int(line[ldx:].split(' ')[0])

                #ldx = line.find(e_search) + len(e_search)
                self.calc_E_z = get_ssv_val(line, e_search, float)#self.max_E_z = float(line[ldx:].split(' ')[0])

                ldx = line.rfind(k_search) + len(k_search)
                self.dev_K = float(line[ldx+1:].split(' ')[0])
            
                tdx = line.find('(', ldx) + 1
                ena_text = 'enabled'
                dis_text = 'disable' # needs to be the same length as enabled, so don't use disabled
                text = line[tdx:tdx + len(ena_text)]
                if text.lower() == ena_text:
                    self.dev_en = True
                elif text.lower() == dis_text:
                    self.dev_en = False
                else:
                    print(text.lower())
                    raise RuntimeError("Parameter line doesn't specify whether devonshire was enabled!")
                found_param_line = True
            elif Emax_search in line:
                # parse 
                # note: the last 
                self.max_E_z = max(get_ssv_val(line, Emax_search, float), self.max_E_z)
            elif "have taken" in line:
                lline = line.lower()
                print(f'found have taken line "{lline}"')
                dur = get_ssv_val(lline, 'have taken', float)
                if lline.startswith('finished matrix elt'):
                    self.mat_elt_dur = dur
                    print(f'found mat_elt_dur {dur}')
                elif lline.startswith('completed stark loop'):
                    self.stk_lop_dur = dur
                    print(f'found stark loop dur {dur}')
        f.close()
        if not found_param_line:
            raise RuntimeError(f"Parameter line not found in {self.log_path}")
        if not math.isfinite(self.max_E_z):
            raise RuntimeError(f"Unable to find maximum Electric field strength in {self.log_path}")
        return self

    def parse_gnd_stark_shift(self, *args, **kwargs):
        fpath = os.path.join(self.path, 'stark_shift_gnd.csv')
        return pd.read_csv(fpath, *args, **kwargs)

    def parse_stark_spect(self, *args, **kwargs):
        fpath = os.path.join(self.path, 'stark_spectrum.csv')
        return pd.read_csv(fpath, *args, **kwargs)

    def get_coeff_dir(self):
        return self.cdir
    
    def has_coeff_dir(self):
        return os.path.exists(self.cdir)

    def parse_state_coeffs(self, E_z, *args, **kwargs):
        csvpath = os.path.join(self.cdir, f'{E_z}.csv')
        return pd.read_csv(csvpath, *args, **kwargs)

    def list_state_csvs(self):
        csvlist = []
        for ent in os.scandir(self.cdir):
            if ent.is_file() and ent.name.lower().endswith('.csv'):
                csvlist.append(ent.path)
        return csvlist
        
    def get_state_ifo_dir(self):
        return self.state_info_dir

    def get_state_ifo_Ez(self, Ez, *args, **kwargs):
        csvpath = os.path.join(self.get_state_ifo_dir(), f'info_Ez_{Ez}.csv')
        return pd.read_csv(csvpath, *args, **kwargs)

    def list_ifo_csvs(self):
        csvlist = []
        for ent in os.scandir(self.get_state_ifo_dir()):
            if ent.is_file() and ent.name.lower().endswith('.csv'):
                csvlist.append(ent.path)
        return csvlist

def find_runs(scandir):
    runlist = []
    for rundir in os.scandir(scandir):
        if matches_run_spec(rundir.name) and rundir.is_dir():
            try:
                runlist.append(aef_run(rundir.path))
            except BaseException as e:
                print(f'directory {rundir} has an invalid run.  See exception for more information')
                print(e)
    return runlist





