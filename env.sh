#!/bin/bash
if [ x$BASH_SOURCE = x -o $BASH_SOURCE == $0 ] ; then
  echo "This script ${0} should be sourced, not run"
  exit 254
fi
scr_src="${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}"
scr_src="$(readlink -f ${scr_src})"
export AEF_DIR="${scr_src%/*}"
echo AeF-hyperfine-structure directory is "${AEF_DIR}"
export PATH="${AEF_DIR}":"$PATH"
pyscr_path=${AEF_DIR}/analysis-plot-scripts
export PYTHONPATH="${pyscr_path}":"${pyscr_path}/src":"$PYTHONPATH"
alias aef_source_env="source $scr_src"

## run plot_state_jbasis.py and log output
function plot_state_jbasis {
  if [ -z $1 -o ! -d $1 -o ! -f $1/out.log ]; then
    echo "specify an aef run directory"
    return
  fi
  olog=/dev/null/invalid_path
  re='^[+-]?[0-9]+([.][0-9]+)?$'
  if ! [[ $2 =~ $re ]] ; then
   olog=$1/state_coeffs/jb_default.log
  else
   olog=$1/state_coeffs/jb_$2.log
  fi
  python3 "$pyscr_path/plot_state_jbasis.py" ${1+"$@"} 2>&1 | tee $olog
}

## run low_state_dumper and log output
function lowstatedump {
  if [ -z $1 -o ! -d $1 -o ! -f $1/out.log ]; then
    echo "specify an aef run directory"
    return
  fi
  out=$1
  $AEF_DIR/low_state_dumper -e 0 -e 500 -e 50000 -l $out/matrix.dat ${1+"$@"}
}
