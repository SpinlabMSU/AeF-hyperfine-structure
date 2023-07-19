
#!/usr/bin/env python3
### plot molecular dipole z-component as a function of E-field


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


#####
arr_Ez = []
arr_Dz = []

if len(sys.argv) < 2:
    print(f"usage: {sys.argv[0]} <dirname>")
    #exit(1)
    #dirname = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-12-162547.3188494\devonshire_info'
    #dirname = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-17-201453.3155610'
    #dirname = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-17-232644.7464130'
    #### ALL 
    dirname = r'C:\Users\nusgart\source\AeF-hyperfine-structure\output\2023-07-19-181153.8494779'
else:
    dirname = sys.argv[1]
    print(dirname)
    print(type(dirname))

top_dir = dirname

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


def doscan(tdir):
    return tdir ## XXX 
    # avoid scan if either 
    date_rgx = '[0-9]+-[0-9]+-[0-9]+'
    number_rgx = '[0-9]+(.[0-9]+)?'
    srgx = fr'.*[/\]{date_rgx}-{number_rgx}([/\](devonshire_info[/\]?)?)?'
    rgx = re.compile(srgx)
    if rgx.match(tdir) != None:
        return tdir
    # otherwise, we're in the directory, need to find latest
    ### TODO actually do this
    return tdir ## do nothing for now

# was devonshire enabled?
def was_dev_en(dir):
    #
    fnam = os.path.join(dirname, "out.log")
    f = open(fnam, 'r')
    for line in f:
        k_search = "K is" # devonshire coupling constant named K, look for enabled
        if "nmax is" in line and k_search in line:
            print(f"Found parameter line \"{line}\"")
            ldx = line.rfind(k_search) + len(k_search)

            dev_K = float(line[ldx+1:].split(' ')[0])
            
            tdx = line.find('(', ldx) + 1
            ena_text = 'enabled'
            dis_text = 'disable' # needs to be the same length as enabled, so don't use disabled
            text = line[tdx:tdx + len(ena_text)]
            if text.lower() == ena_text:
                return True, dev_K
            elif text.lower() == dis_text:
                return False, dev_K
            else:
                raise RuntimeError("Parameter line doesn't specify whether devonshire was enabled!")
    raise RuntimeError(f"Devonshire status not found in {fnam}")
dev_en, dev_K = was_dev_en(dirname)
# scan for devonshire_info and use that 
tdir = os.path.join(dirname, "devonshire_info")
if os.path.exists(tdir) and os.path.isdir(tdir):
    dirname = tdir
else:
    dirname = doscan(dirname)

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
    arr_Ez.append(Ez)
    arr_Dz.append(dz)

print(arr_Ez)
print(arr_Dz)

Ez = np.array(arr_Ez)
Dz = np.array(arr_Dz)

fig = plt.figure(figsize=(13.66, 7.68))
status_txt = f"enabled, K={dev_K}" if dev_en else "disabled"
run = os.path.split(top_dir)[1]
title_text = f"Molecular Dipole-Axis Z vs Externally applied Electric Field (with Devonshire {status_txt})\nrun {run}"
plt.title(title_text)
plt.xlabel('Externally applied electric field (V/cm)')
plt.ylabel('Molecular Dipole-Axis z component (unitless)')
plt.plot(Ez, Dz, 'ro')
plt.savefig(os.path.join(top_dir, 'mdz_stark_plot.png'))
plt.show()
