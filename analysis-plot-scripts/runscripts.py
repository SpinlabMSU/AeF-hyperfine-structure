#!/usr/bin/env python3
## runscripts.py -- . 
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
import sys
import os.path
import subprocess

exec_path = sys.executable
script_dir = os.path.split(os.path.realpath(__file__))[0]
base_path = os.path.split(script_dir)[0]

outdir = os.path.join(base_path, "output")

rundir = None
rundir2 = None

if len(sys.argv > 1):
    rundir = sys.argv[1]
elif len(sys.argv > 2):
    rundir2 = sys.argv[2]

if rundir == None:
    pass

def run_script(sc_nam, *args):
    script = os.path.join(base_path, sc_nam + ".py")
    return subprocess.run([exec_path, script, *args])



run_script()