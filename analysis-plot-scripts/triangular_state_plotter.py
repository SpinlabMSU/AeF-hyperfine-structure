#!/usr/bin/env python3
## triangular_state_plotter.py -- implements code to plot the j-basis.
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

def draw_pan(ax, x_start, x_end, y, h, txt, *args, **kwargs):
    ax.hlines(y, x_start, x_end, *args, **kwargs)
    ax.vlines([x_start, x_end], y, y+h, **kwargs)
    x_txt = (x_start + x_end) / 2
    ax.annotate(txt, (x_txt, y + 0.1), ha='center')

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