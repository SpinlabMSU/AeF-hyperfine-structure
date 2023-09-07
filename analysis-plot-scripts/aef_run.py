## aef_run.py: 
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

class aef_run(object):
    
    def __init__(self, dir_path):
        self.path = dir_path
        self.run = os.path.split(dir_path)[1]
        self.timestamp = datetime.strptime(self.run, "%Y-%m-%d-%H%M%.%f")

        self.log_path = os.path.join(self.path, 'out.log')

        self.valid = self.check_valid()
        self.dev_en = False
        self.dev_K = math.nan
        self.nmax = -1 # n is the orbitorotational quantum number
        self.max_E_z = math.nan

        if self.valid:
            self.parse_log()

    def check_valid(self):
        run_exists = os.path.exists(self.path) and os.path.isdir(self.path)
        log_exists = os.path.exists(self.log_path) and os.path.isfile(self.log_path)
        self.valid = run_exists and log_exists
        return self.valid

    def parse_log(self):
        f = open(self.log_path, 'r')
        for line in f:
            n_search = "nmax is "
            e_search = "E_z is " # externally applied
            k_search = "K is" # devonshire coupling constant named K, look for enabled
            if "nmax is" in line and k_search in line:
                print(f"Found parameter line \"{line}\"")

                ldx = line.find(n_search) + len(n_search)
                self.nmax = int(line[ldx:].split(' ')[0])

                ldx = line.find(e_search) + len(e_search)
                self.max_E_z = float(line[ldx:].split(' ')[0])

                ldx = line.rfind(k_search) + len(k_search)
                dev_K = float(line[ldx+1:].split(' ')[0])
            
                tdx = line.find('(', ldx) + 1
                ena_text = 'enabled'
                dis_text = 'disable' # needs to be the same length as enabled, so don't use disabled
                text = line[tdx:tdx + len(ena_text)]
                if text.lower() == ena_text:
                    self.dev_en = True
                elif text.lower() == dis_text:
                    self.dev_en = False
                else:
                    raise RuntimeError("Parameter line doesn't specify whether devonshire was enabled!")
                f.close()
                break
            #elif line.
        raise RuntimeError(f"Devonshire status not found in {self.log_path}")

    def parse_gnd_stark_shift(self, *args, **kwargs):
        fpath = os.path.join(self.path, 'stark_shift_gnd.csv')
        return pd.read_csv(fpath, *args, **kwargs)

    def parse_stark_spect(self, *args, **kwargs):
        fpath = os.path.join(self.path, 'stark_spectrum.csv')
        return pd.read_csv(fpath, *args, **kwargs)

def find_runs(scandir):
    runlist = []
    for rundir in os.scandir(scandir):
        if matches_run_spec(rundir.name) and rundir.is_dir():
            runlist.append(aef_run(rundir.path))
    return runlist





