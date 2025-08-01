# Compiling AeF-hyperfine-structure

This document provides instructions on how to compile the AeF-hyperfine-structure toolkit.  Please note that no
installation is required -- you can directly run the executables immediately after compiling them.

Specific instructions will be provided for Linux and Windows, but in general .  Please note that Mac OS is not a
supported platform at the present time due to lacking a macbook to test on and also as a result of CUDA not being
supported on macOS.

## Dependencies
* A C++ 23 compliant compiler (G++ 13 and MSVC v19.37 / VS 2022 version 17.7 implement enough of C++ 23 to work)
* Eigen v3.4.0 -- included
* gcem -- included
* CUDA Toolkit -- specifically cuBLAS and cuSOLVER
* CERN ROOT v6 -- only tested with version 6.32, but any recent version of ROOT 6 should work. Please note that ROOT 7,
when it is released, may not be backwards compatible.
* zlib
* GNU GSL
* (Optionally) Intel MKL

## Compiling on Windows
Install all of the required dependencies and use the visual studio project.  You should install visual studio 2022
(or a later version, though this is untested) from https://visualstudio.microsoft.com/downloads/.  The free community
option suffices. Make sure to install the "Desktop development with C++" option. To install CUDA, follow NVIDIA's
instructions at https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/. It is suggested to use vcpkg
to install the other dependencies.

## Compiling on Linux

Install all of the required dependencies.  On debian-based distributions, most dependencies can be installed by
```
sudo apt install build-essential libgsl-dev zlib1g-dev
```
You will then need to install NVIDIA's CUDA toolkit following the instructions at https://docs.nvidia.com/cuda/cuda-installation-guide-linux/