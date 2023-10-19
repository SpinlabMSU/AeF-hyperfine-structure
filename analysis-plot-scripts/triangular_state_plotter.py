### Triangular_state_plotter
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
import baf_state

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
    bidx = 4 * n * n
    num_states = 4 * n + 8
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
#def get_indicies_njf(n, j, f):
#    return ur mom

class nlevel:
    def __init__(self, n):
        assert(isinstance(n,int))
        self.n = n
        self.bidx = 4 * n * n
        self.num_states = 8 * n + 4
        self.indicies = np.arange(self.bidx, self.bidx + self.num_states, 1, dtype=np.int64)
        self.jfs = []
        if n > 0:
            jm = n - 0.5
            self.jfs += ((jm, jm - 0.5), (jm, jm + 0.5))
        jp = n + 0.5
        self.jfs += ((jp, jp - 0.5), (jp, jp + 0.5))
        self.jfs = tuple(self.jfs)
        self.gidx_jfs = tuple(range(4)) 
        self.jf_n_mfs = tuple(2*f+1 for (j, f) in self.jfs)
        self.max_mfs = 2 * n + 3 ## is 2 * (n + 0.5 + 0.5) + 1 == 2 * (n + 1) + 1

    def gidx(self, j, f):
        """
        Returns the "group-index" of the jf pair in this n level.
        """
        for idx in range(len(self.jfs)):
            (jj, ff) = self.jfs[idx]
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
        height = (m_f + f) + self.base_offset(gidx)
        return (gidx, height)

    def draw_boxes(self, ax, ofs=0):
        for idx in self.indicies:
            st = baf_state.state_from_index(int(idx))
            x, y = self.get_box_pos(st.j, st.f, st.m_f)
            x += ofs
            print(st, x, y)
            ax.plot(x, y, 'ro')
            ax.annotate(f'|{st.n},{st.j}\n{st.f},{st.m_f}>', (x, y))


if __name__ == '__main__':
    print('running')
    fig = plt.figure()
    njs = 0
    for n in range(4):
        nl = nlevel(n)
        nl.draw_boxes(plt.gca(), njs)
        prev_njs = njs
        njs += 4 if (n > 0) else 2
        plt.hlines(-2, prev_njs, njs - 0.5, color='k', linestyle='-', label=f'n={n}')
        dx = ((njs - 0.5) - (prev_njs)) / 2
        x = prev_njs + 0.5
        plt.text(x + dx, -1.5, f'n={n}', ha='right', va='center')
    plt.show()