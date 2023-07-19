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