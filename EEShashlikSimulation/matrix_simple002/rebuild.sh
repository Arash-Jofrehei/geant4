echo "Remember to set up the environment first with setup_geant4.sh"
if [ $HOSTNAME = "t3ui17" ]; then
    echo "Rebuilding on PSI"
    cd /shome/tklijnsm/testbeam/geant4/EEShashlikSimulation/H4OpticalSmall
    rm -rf cmake
    mkdir cmake
    cd cmake
    cp ../rebuild.sh rebuild.sh
    cp ../recompile.sh recompile.sh
    cp ../runThomas.mac runThomas.mac
    cmake -DGeant4_DIR=/shome/tklijnsm/geant4/geant4.10.01.p02-install/lib64/Geant4-10.1.2/ ..
    make
else
    echo "Rebuilding on lxplus"
    cd /afs/cern.ch/user/t/tklijnsm/Work/geant4/EEShashlikSimulation/H4OpticalSmall
    rm -rf cmake
    mkdir cmake
    cd cmake
    cp ../rebuild.sh rebuild.sh
    cp ../recompile.sh recompile.sh
    cp ../runThomas.mac runThomas.mac
    cmake ..
    make
fi
