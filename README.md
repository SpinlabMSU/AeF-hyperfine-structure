# AeF-hyperfine-structure: Alkaline-Earth Monofluoride molecular hyperfine structure calculator
AeF-hyperfine-strucutre is a toolkit for computing the rotational-hyperfine structure of alkaline-earth monofluoride
molecules in their $^2\Sigma+$ electronic-vibrational ground states.  Instructions on how to compile the  software are
provided in `COMPILING.MD`, and the license (GPLv3) is provided in `COPYING`.  The purpose of this toolkit is
specifically calculating spectra and matrix elements useful for determining measurement schemes for CP-violation searches
using alkaline-earth monofluorides in solid noble gas matricies, especially for nuclear schiff moment searches using
radium-225 monofluoride embedded in a solid argon matrix.  Towards this


## Provided programs
* AeF-hyperfine-structure: this program calculates the rotational-hyperfine spectrum of an alkaline-earth monofluoride 
across a range of electric field strengths either in the gas phase or embedded in a solid medium.  Currently, only 
$^{138}\mathrm{BaF}$ is supported, but work is ongoing to enable calculations on $^{225}\mathrm{RaF}$ as well.
* LowStateDumper -- this program outputs selected expectation values and matrix elements of the lowest set of energy eigenstates
* PerturbationAnalyzer: this program performs first-order perturbative calculations of a selected set of "interaction" operators
* operator\_visualizer: this program dumps
* StarkDiagonalizer -- this program performs similar spectrum calculations to AeF-hyperfine-structure but only includes the Stark
interaction in the Hamiltonian.  Its output is useful for comparison purposes.
* GenerateHamiltonianFiles -- deprecated, do not use
* NoStark\_HyperfineTester -- this is a "playground" testing program that is only useful for debugging SpinlabHyperfineLib

## Using the AeF-hyperfine-structure toolkit
The flow



## 

## A Brief theoretical overview of the 

$\newcommand{\ket}[1]{\left|{#1}\right\rangle}$
$\newcommand{\bra}[1]{\left\langle{#1}\right|}$
$\newcommand{\abs}[1]{\left\vert{#1}\right\vert}$
$\newcommand{\braket}[2]{\left\langle{#1}\middle|{#2}\right\rangle}$

This code here is a testbed for computing the hyperfine structure of $^{138}BaF$ molecules in their $^2\Sigma$ electronic and vibrational ground states.  This is based on PRA 98, 032513 (2018) using 

Under these circumstances, the state of each molecule can be described using three coupled angular momenta:
* $\vec{I}$: total nuclear spin ($I = \frac{1}{2}$ always since $^{138}Ba$ has $I=\frac{1}{2}$ and $^{19}F$ has $I=0$)
* $\vec{S}$: total electron spin ($S = \frac{1}{2}$ in the electronic ground state)
* $\vec{N}$: molecular rotational angular momentum ($n\in\mathbb{Z}$)

The total angular momentum of the molecule is denoted $\vec{F}=\vec{I}+\vec{S}+\vec{N}$ with quantum numbers $f$ and $m_f$.  Since there are three angular momenta, there are several possible bases

## J-basis
One possible basis couples $S$ and $N$ first to form $\vec{J} = \vec{N} + \vec{S}$ ($j=n\pm\frac{1}{2}$), and then couples I to $J$ to form $F$.  This will be used as the "default" basis since $S$ and $N$ couple more strongly to each other than to $I$.

For a given $n\neq0$, the quantum number $j$ has two possible values $j^+(n) = n + \frac{1}{2}$ and $j^-(n) = n - \frac{1}{2}$.  For $n=0$, only $j^+(0)=\frac{1}{2}$ is defined.

Similarly, $f$ has two possible values for a given $j$, $f^+(j) = j + \frac{1}{2}$ and $f^-(j) = j - \frac{1}{2}$.  Note that $j$ is always a half-integer, so both $f^\pm$ are defined for all $j$. 

Thus, in the $J$-basis a molecular state can be described as $\ket{i(sn)jfm_f}$.  This can be abbreviated to $\ket{njfm_f}$ since $i$ and $s$ are always $\frac{1}{2}$.

## G-basis
Under certain circumstances, it is useful to couple $I$ and $S$ first instead of $S$ and $N$.  In these circumstances, a new angular momentum $\vec{G} = \vec{I} + \vec{S}$ is used instead.  It has quantum number $g = 0, 1$ with the usual singlet-triplet set.

These two bases are related by the wigner 6j symbol --
$$ \braket{i(sn)jf}{(is)gnf} = \xi'()$$

## 3M-basis
There is also an uncoupled basis $\ket{im_i,sm_s,nm_n}$.  Since $f$ is not diagonal in this basis, this is only really useful in evaluating the "tensor" part of the hyperfine Hamiltonian.


## Hamiltonian in free space
The effective Hamiltonian in vacuum (with possible electric fields can be described as the sum of three parts: a rotational Hamiltonian, a Stark shift, and a hyperfine shift:
$$ H = H_{rot} + H_{st} + H_{hfs} $$

### Rotational Hamiltonian
The rotational Hamiltonian is $$ H_{rot} = BN^2 - DN^4 + \gamma\vec{N}\cdot\vec{S} + \delta N^2 \vec{N}\cdot\vec{S}$$  Note that this is diagonal in the $J-basis$, and 


## J-basis index
It is often useful to assign a natural number index to each element of the most frequently used basis.  For example, this makes it easy to efficiently represent operators as n-d matricies.


## References
1. PRA 98, 032513 (2018) (EDM3 proposal paper)
2. J. Chem. Phys. 105, 7412 (1996).
3. 

## Licensing information
AeF-hyperfine-structure is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

AeF-hyperfine-structure is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
AeF-hyperfine-structure. If not, see <https://www.gnu.org/licenses/>.
