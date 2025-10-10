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
from sympy import KroneckerDelta
import fractions
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
class HConstants:
    ## sourced from PRA 98, 032513 (2018)
    B:ClassVar[rational] = 6743.9586 # MHz
    D:ClassVar[rational] = 5.5296 * 1e-3 # kHz --> MHz
    gamma:ClassVar[rational] = 80.954
    delta:ClassVar[rational] = 0.111 *1e-3 # kHz --> MHz

    b:ClassVar[rational] = 63.509 # MHz
    c:ClassVar[rational] =  8.224 # MHz 
    mu_e: ClassVar[rational] = 3.170 #D


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
        ndx = 4*self.n*self.n
        jdx = 4*self.n if self.j > self.n else 0
        f1dx = 2*self.j if self.f > self.j else 0
        fdx = 0
        mdx = self.m_f + self.f
        idx = ndx + jdx + fdx + mdx
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
    
    def Hrot(self, B=HConstants.B, D=HConstants.D,
             gamma=HConstants.gamma, delta=HConstants.delta):
        nsq = self.Nmag()
        jsq = self.Jmag()
        ssq = self.Smag()
        nds = self.n_dot_s()
        return B*nsq - D * nsq * nsq + (gamma + delta * nsq) * nds
    
    def tostr(self):
        ## intentionally not overrriding __repr__ or __str__
        return f"|n={self.n},j={self.j},f={self.f},m_f={self.m_f}>"

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
    n = isqrt(rem >> 2)
    rem -= 4 * n * n
    j = n + 0.5 if rem >= 4 * n else n - 0.5
    if rem >= 4 * n:
        rem -= 4 * n
    f = j + 0.5 if rem >= 2 * j else j - 0.5
    if rem >= 2 * j:
        rem -= 2 * j
    m_f = rem - f
    return HyperfineState(n, j, f, m_f)

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