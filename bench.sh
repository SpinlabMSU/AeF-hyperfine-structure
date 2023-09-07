#!/bin/bash
# benchmark performance for nmax ranging from 4 to 40

for nmax in `seq 4 40`; do
  ./x64/Rnodev/AeF-hyperfine-structure.exe -n $nmax --print_extras=false
done

for nmax in `seq 4 40`; do
  ./x64/Rdeven/AeF-hyperfine-structure.exe -n $nmax --print_extras=false
done
