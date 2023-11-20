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


def angular_H_st(n,j,f,mf, nn,jj,ff,mff):
    """
    This code is probably broken and is included only 
    """
    ## only includes the
    ## \vec{F} == total angular momentum
    # f --> angular state
    ## \vec{I} == total nuclear spin
    # I think this formula assumes i == 1/2 
    ## \vec{S} == total electron spin
    # 
    ## \vec{J} == \vec{F} - \vec{I}
    # j describes ---> should be
    ## \vec{N} == \vec{J} - \vec{S} : molecular angular momentum
    # n describes the "rotational state" --> should be integer 
    if (mf != mff):
        return 0
    xi_factors = xi(f,ff) * xi(j,jj) * xi(n,nn) * xi(mf, mff)
    a3j_s = wigner_3j(f, 1, ff, -mf, 0, mf) * wigner_3j(n,1,nn,0,0,0)
    
    if a3j_s == 0:
        return 0
    
    a6j_s = wigner_6j(f, 1, ff, jj, half, j)*wigner_6j(j,1,jj,nn,half,n)
    kds = kron_delt(mf, mff)
    phase = sympy.Pow(-1, 1-mf)
    return phase*kds * xi_factors * a3j_s * a6j_s



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
    f: rational
    m_f: rational
        
        
    
    def index(self):
        ndx = 4*self.n*self.n
        jdx = 4*self.n if self.j > self.n else 0
        fdx = 2*self.j if self.f > self.j else 0
        mdx = self.m_f + self.f
        idx = ndx + jdx + fdx + mdx
        assert(idx == int(idx))
        return int(idx)
    
    def Nmag(self):
        return self.n * (self.n + 1)
    
    def Jmag(self):
        return self.j * (self.j + 1)
    
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
    
    
    def _H_hfs_scalar(self, other, b):
        ## scalar part of the molecular hyperfine structure
        if (other.n != self.n or other.f != self.f or other.m_f != self.m_f):
            return 0
        n = self.n; f = self.f; j = self.j; jp = other.j
        
        ## 3 / 4 * b * xi_{jj'}
        coeff = 3 * half * half * b * xi(j,jp)
        hg0 = w6j(half,half,0, n,f,jp) * w6j(half,half,0, n,f,j)
        hg1 = w6j(half,half,1, n,f,jp) * w6j(half,half,1, n,f,j)
        
        return coeff * (hg1 - hg0)
    
    def _H_hfs_tensor(self, other, c):
        j = self.j; jp = other.j; f = self.f; fp = other.f
        m_f = self.m_f; m_fp = other.m_f
        n = self.n
        
        if self.n != other.n or self.f != other.f or self.m_f != other.m_f:
            return 0
        
        htensor = 0
        prf = c * abs(xi(f, fp)) * xi(j, jp) * parity(-2*m_f)
        i = s = half
        for m_j in np.arange(-j, j+1, 1):
            for m_jp in np.arange(-jp, jp+1, 1):
                for m_i in (-half, half):
                    for m_s in (-half, half):
                        for m_n in range(-n, n+1, 1):
                            ##
                            w0 = w3j(i,jp,fp,m_i,m_jp,-m_fp)
                            w1 = w3j(s,n,jp,m_s,m_n,-m_jp)
                            w2 = w3j(i,j,f,m_i,m_j,-m_f)
                            w3 = w3j(s,n,j,m_s,m_n,-m_j)
                            htensor += w0*w1*w2*w3*m_i*m_s*parity(-m_j-m_jp)
        return prf * htensor
    
    def Hhfs(self, other, b=HConstants.b, c=HConstants.c):
        # b * \vec{I}\cdot\vec{S} + cI_{z}S_{z}
        ## need to sum over all g
        ## scalar part:
        hscalar = self._H_hfs_scalar(other, b)
        # tensor part:
        htensor = self._H_hfs_tensor(other, c)
        return hscalar + htensor
    
    def Hst(self, other, E_z, mu_e=HConstants.mu_e):
        coeff = angular_H_st( self.n,  self.j,  self.f,  self.m_f,
                             other.n, other.j, other.f, other.m_f)
        return mu_e * E_z * coeff
    
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
                    state = HyperfineState(n,j,f,m_f)
                    states.append(state)
                    m_f += 1
    return states

if __name__ == '__main__':
    nmax = 8
    hsts = make_hyperfine_states(nmax)
    pprint.pprint((hsts))
    print(f'have {len(hsts)} hyperfine states with nmax={nmax}')