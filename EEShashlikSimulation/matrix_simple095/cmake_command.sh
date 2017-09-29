mkdir cmake
source ./setup_g4.sh
source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/g++
export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/bin/gcc
source /afs/cern.ch/sw/lcg/external/geant4/10.1.p01/x86_64-slc6-gcc49-opt/CMake-setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
echo "Environment set up for lxplus"
cd cmake
mkdir ntuples
#cmake -DGeant4_DIR=/Users/pandolf/local/geant4.10.00.p01-install/lib/Geant4-10.0.1/ ..
cmake -DGeant4_DIR=/afs/cern.ch/sw/lcg/external/geant4/10.1.p01/x86_64-slc6-gcc49-opt/lib64/Geant4-10.1.1/ ..
make -j 8
./runEEShashlik -m ../../run2.mac