#!/bin/bash
cd cmake/
if [[ $HOSTNAME = "t3ui19" ]]; then
    # Set paths to the extra Geant4 data; works by default on lxplus
    export G4LEDATA=/shome/tklijnsm/geant4/geant4.10.01.p02-install/share/Geant4-10.1.2/data/G4EMLOW6.40
    export G4LEVELGAMMADATA=/shome/tklijnsm/geant4/geant4.10.01.p02-install/share/Geant4-10.1.2/data/PhotonEvaporation3.1
    export G4EMLOW=/shome/tklijnsm/geant4/geant4.10.01.p02-install/share/Geant4-10.1.2/data/G4EMLOW6.40
    export G4REALSURFACEDATA=/shome/tklijnsm/geant4/geant4.10.01.p02-install/share/Geant4-10.1.2/data/RealSurface1.0

    # TEST IF THIS IS ACTUALLY NEEDED - outdated cmake seems to work okay actually
    #export PATH=/afs/psi.ch/sys/psi.x86_64_slp6/Programming/cmake/2.8.12.2/bin:$PATH
    #which cmake /afs/psi.ch/sys/psi.x86_64_slp6/Programming/cmake/2.8.12.2/bin/cmake


    source /afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/setup.sh
    source /afs/cern.ch/sw/lcg/external/geant4/10.0.p01/x86_64-slc6-gcc47-opt/CMake-setup.sh
    export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/g++
    export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/gcc



    echo "Environment set up for PSI"
    echo "Remember to also set up CMSSW"
else
#    source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.32/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh    
    source /afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/setup.sh
    source /afs/cern.ch/sw/lcg/external/geant4/10.0.p01/x86_64-slc6-gcc47-opt/CMake-setup.sh
    export CXX=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/g++
    export CC=/afs/cern.ch/sw/lcg/contrib/gcc/4.7/x86_64-slc6/bin/gcc
    echo "Environment set up for lxplus"
fi
