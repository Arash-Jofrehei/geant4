#!/bin/bash

### SOURCE ROOT
echo "Sourcing ROOT v5.34.18 py26"
#source /swshare/ROOT/root_v5.34.18_slc5_amd64_py26_pythia6/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/external/geant4/10.0.p01/x86_64-slc6-gcc47-opt/CMake-setup.sh
export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/g++
export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/gcc
source /mnt/t3nfs01/data01/swshare/ROOT/root_v5.34.18_slc6_amd64_py26_pythia6/bin/thisroot.sh
#source /swshare/ROOT/root_v5.34.32_precompiled/root/bin/thisroot.sh 


