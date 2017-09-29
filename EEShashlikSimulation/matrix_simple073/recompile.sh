if [ $HOSTNAME = "t3ui17" ]; then
    echo "Recompiling on PSI"
    cd /shome/tklijnsm/testbeam/geant4/EEShashlikSimulation/H4OpticalSmall/cmake
    cmake -DGeant4_DIR=/shome/tklijnsm/geant4/geant4.10.01.p02-install/lib64/Geant4-10.1.2/ ..
else
    echo "Recompiling on lxplus"
    cd /afs/cern.ch/user/t/tklijnsm/Work/geant4/EEShashlikSimulation/H4OpticalSmall/cmake
    cmake ..
fi
make
