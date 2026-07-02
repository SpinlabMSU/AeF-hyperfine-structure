#!/usr/bin/env python3
## raf_state.py -- implements the jf-basis for tracking the rotational-hyperfine
# structure of 225RaF
## This does not implement any matrix elements, unlike baf_state.py since they're
# broken and unusable anyways
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
## This code implements
import numpy as np
import matplotlib.pyplot as plt
import sympy as sympy
import math
import numba
from sympy.physics.wigner import *
from fractions import Fraction
from dataclasses import dataclass
import pprint
import sys
from typing import ClassVar
#import warnings

use_sympy = not True

if use_sympy:
    def make_rational(a, b):
        return Fraction(a, b)
else:
    def make_rational(a, b):
        return a / b
half = make_rational(1, 2)
rational = Fraction if use_sympy else float



def kron_delt(i,j):
    return (i == j)+0

def parity_py(x):
    return (-1)**x

def parity_sympy(x):
    return sympy.Pow(-1, x)


parity = parity_sympy if use_sympy else parity_py
def xi_sympy(s, ss):
    s = sympy.S(s)
    ss = sympy.S(ss)
    n_s = 2 * s + 1
    n_sprime = 2 * ss + 1
    return sympy.Pow(sympy.S(-1), s+ss) * sympy.sqrt(n_s * n_sprime)

def xi_np(s, sprime):
    n_s = 2 * s + 1
    n_sprime = 2 * sprime + 1
    return (-1)**(s+sprime) * np.sqrt(n_s * n_sprime)

@numba.njit
def xi_fast(s, sprime):
    n_s = 2 * s + 1
    n_sprime = 2 * sprime + 1
    return (-1+0j)**(s+sprime) * math.sqrt(n_s * n_sprime)

xi = xi_sympy if use_sympy else xi_fast
w3j = wigner_3j
w6j = wigner_6j

@dataclass
class HyperfineState:
    # Quantum operators are F,I,S,J,N
    # quantization axis is along the Ba-F 
    # \vec{F} = total angular momentum
    # \vec{I} = total nuclear spin, i = 1/2
    # \vec{S} = total electron spin, s = 1/2
    # \vec{J} = \vec{F} - \vec{I} = 
    n: rational  ## all in units of hbar, of course
    j: rational
    f_1:rational
    f: rational
    m_f: rational
    
    def index(self):
        ndx = 8*self.n*self.n
        jdx = 8*self.n if self.j > self.n else 0
        f1dx = 4*self.f1 if self.f > self.j else 0
        fdx = 2*self.f1 if self.f > self.f1 else 0
        mdx = self.m_f + self.f
        idx = ndx + jdx +f1dx + fdx + mdx
        assert(idx == int(idx))
        return int(idx)
    
    def Nmag(self):
        return self.n * (self.n + 1)
    
    def Jmag(self):
        return self.j * (self.j + 1)

    def F1mag(self):
        return self.f_1 * (self.f_1 + 1)
    
    def Fmag(self):
        return self.f * (self.f + 1)
    
    def Smag(self):
        return half * (half + 1)
    
    def n_dot_s(self):
        return half * (self.Jmag() - self.Nmag() - self.Smag())
    
    def tostr(self):
        ## intentionally not overrriding __repr__ or __str__
        return f"|n={self.n},j={self.j},f1={self.f1},f={self.f},m_f={self.m_f}>"

def isqrt(n):
    if n > 0:
        x = 1 << (n.bit_length() + 1 >> 1)
        while True:
            y = (x + n // x) >> 1
            if y >= x:
                return x
            x = y
    elif n == 0:
        return 0
    else:
        raise ValueError("square root not defined for negative numbers")

def state_from_index(idx):
    rem = idx
    n = isqrt(rem >> 3)
    rem -= 8 * n * n
    j = n + half if rem >= 8 * n else n - half
    if rem >= 8 * n:
        rem -= 8 * n
    f1 = j + half if rem >= 4 * j else j - half
    if rem >= 4 * j:
        rem -= 4 * j
    f = f1 + half if rem >= 2 * f1 else f1 - half
    if rem >= 2 * f1:
        rem -= 2 * f1
    m_f = rem - f
    return HyperfineState(n, j, f1, f, m_f)

def make_hyperfine_states(nmax):
    states = []
    for n in range(nmax+1):
        for j in [n-half, n+half]:
            if j < 0: continue
            for f in [j-half,j+half]:
                if f < 0: continue
                m_f = -f
                while m_f <= f:
                    state = HyperfineState(n,j,f_1,f,m_f)
                    states.append(state)
                    m_f += 1
    return states

"""
For n=0 should have a plot roughly looking like
   
   xxx                         --------------------
   xxx                         | j=0.5,f=1,m_f=+1 |
   xxx                         --------------------
   xxx   --------------------  --------------------
   xxx   | j=0.5,f=1,m_f=+1 |  | j=0.5,f=1,m_f= 0 |
   xxx   --------------------  --------------------
   xxx                         --------------------
   xxx                         | j=0.5,f=1,m_f=-1 |
   xxx                         --------------------
       |     j = 0.5,f = 0        j=0.5, f = 1       |
       |---------------------------------------------|
                         n = 0
Where the states
"""


## general design: 


def get_indicies_for_n(n):
    bidx = 8 * n * n
    num_states = 8 * n + 16
    return np.arange(bidx, bidx + num_states, 1)

def jfs_for_n(n):
    jfs = []
    # j- for n > 0
    if n > 0:
        jm = n - 0.5
        jfs += ((jm, jm - 0.5), (jm, jm + 0.5))
    # j+ for all n
    jp = n + 0.5
    jfs += ((jp, jp - 0.5), (jp, jp + 0.5))
    return jfs

def draw_pan(ax, x_start, x_end, y, h, txt, *args, **kwargs):
    ax.hlines(y, x_start, x_end, *args, **kwargs)
    ax.vlines([x_start, x_end], y, y+h, **kwargs)
    x_txt = (x_start + x_end) / 2
    ax.annotate(txt, (x_txt, y + 0.1), ha='center')

class nlevel:
    def __init__(self, n):
        assert(isinstance(n,int))
        self.n = n
        self.bidx = 8 * n * n
        self.num_states = 16 * n + 8
        self.indicies = np.arange(self.bidx, self.bidx + self.num_states, 1, dtype=np.int64)
        self.jffs = []
        if n > 0:
            jm = n - 0.5
            if (jm > 0):
                f1m = jm - 0.5
                self.jffs += ((jm, f1m, f1m - 0.5), (jm, jm + 0.5))
            self.jffs += ((jm, jm - 0.5), (jm, jm + 0.5))
        jp = n + 0.5
        self.jffs += ((jp, jp - 0.5), (jp, jp + 0.5))
        self.jffs = tuple(self.jffs)
        self.gidx_jfs = tuple(range(4)) 
        self.jf_n_mfs = tuple(2*f+1 for (j, f) in self.jffs)
        self.max_mfs = 2 * n + 3 ## is 2 * (n + 0.5 + 0.5) + 1 == 2 * (n + 1) + 1

    def gidx(self, j, f):
        """
        Returns the "group-index" of the jf pair in this n level.
        """
        for idx in range(len(self.jffs)):
            (jj, ff) = self.jffs[idx]
            if j == jj and f == ff:
                return idx
        return -1
    def base_offset(self, gidx):
        # This would always be an integer even without the // since f is always an integer (for any such system of
        # coupled angular momenta, all values of f must be the same "half-parity")
        offs = (self.max_mfs - self.jf_n_mfs[gidx]) // 2
        return offs

    def get_box_pos(self, j, f, m_f):
        gidx = self.gidx(j, f)
        height = m_f #(m_f + f) + self.base_offset(gidx)
        return (gidx, height)

    def draw_boxes(self, ax:plt.Axes, ofs:float=0, func=None, no_text=False):
        """
        Misnomer -- draws the states associated with this n level, including the "pans"
        """
        x = 0
        y = 0
        for idx in self.indicies:
            st = baf_state.state_from_index(int(idx))
            x, y = self.get_box_pos(st.j, st.f, st.m_f)
            x += ofs
            print(st, x, y)
            if func != None:
                func(ax, st, idx, x, y)
            else: ax.plot(x, y, 'ro')
            #ax.annotate(f'|{st.n},{st.j}\n{st.f},{st.m_f}>', (x, y))
            if not no_text:
                ax.annotate(f'|{st.f},{st.m_f}>', (x, y + 0.1), ha='center')
            if st.m_f == -st.f:
                draw_pan(ax, x - 0.45, x + 0.45, y - 0.5, 2*st.f+1, f'f={st.j}', color='m', linestyle='-')
                ax.annotate(f'f={st.f}', (x, y + 2 * st.f + 1), ha='center')
                #ax.annotate(f'f={st.f}', (x, y - 0.5), ha='center')
                #plt.hlines(y - 0.6, x - 0.4, x + 0.4, color='m')
                if st.f > st.j:
                    # j+
                    draw_pan(ax, x - 1.5, x + 0.5, y - 0.9, 0.4, f'j={st.j}', color='b', linestyle='-')
                    #plt.hlines(y - 0.9, x - 1.25, x + 0.25, color='b', linestyle='-')
                    #ax.annotate(f'j={st.j}', (x - 0.5, y - 0.8), ha='center')
        x_start = ofs - 0.8; x_end = x + 0.8
        y = -y
        draw_pan(ax, x_start, x_end, y - 1.2, 1, f'n={self.n}', color='g', linestyle='-')
        x_txt = ofs + (x - ofs) / 2.0
        #ax.annotate(f'n={self.n}', (x_txt, y - 1.1), ha='center')
        #plt.hlines(y - 1.2, x_start, x_end, color='g', linestyle='-')

if __name__ == '__main__':
    print('running')
    fig = plt.figure(figsize=(19.2, 16.8))
    njs = 0
    for n in range(4):
        nl = nlevel(n)
        nl.draw_boxes(plt.gca(), njs)
        prev_njs = njs
        njs += 5 if (n > 0) else 3
        #plt.hlines(-5, prev_njs, njs - 0.5, color='k', linestyle='-', label=f'n={n}')
        dx = ((njs - 0.5) - (prev_njs)) / 2
        x = prev_njs + 0.5
        #plt.text(x + dx, -4.5, f'n={n}', ha='right', va='center')
    plt.title('J-basis basis-state plot')
    plt.ylabel('m_f')
    plt.xlabel('Arbitrary')
    plt.savefig('j-basis-plot.png')
    plt.show()


if __name__ == '__main__':
    nmax = 40
    hsts = make_hyperfine_states(nmax)

    # Tested good up til nmax = 40
    ## should be good for all nmax, but we only need up to 40
    for idx in range(len(hsts)):
        hst = hsts[idx]
        assert(idx == hst.index())
        assert(state_from_index(idx).index() == idx)
    pprint.pprint((hsts))
    print(f'have {len(hsts)} hyperfine states with nmax={nmax}')