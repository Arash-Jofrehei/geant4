#! /usr/bin/env python
import os
import sys
import time
import re
import math
from optparse import OptionParser


def get_string(x):
  thestring = str(int(x))
  if ( math.modf(x)[0] > 0. ) :
    thestring = str(int(math.modf(x)[1]) ) + "p" + str(int(10.*(math.modf(x)[0])) )
  return thestring






parser = OptionParser()
parser.add_option("-q","--queue",dest="queue",default="8nh")
parser.add_option("-l","--label",dest="batchName",default="test")
parser.add_option("-n","--nLayers",dest="nLayers",default=20)
parser.add_option("-a","--act",dest="activeLayerThickness",default=5)
parser.add_option("","--abs",dest="absorberLayerThickness",default=10)
parser.add_option("","--transv",dest="transverseCellSize", default=20)
parser.add_option("-i","--impact",dest="impactPosition", default=0)
parser.add_option("-e","--energy",dest="energy", default=0.75)

(options,args)=parser.parse_args()


#if options.batchName=="nonamenoname" :
#  print( "please provide a label with the -l or --label option" )
#  sys.exit(1)




dir = "/afs/cern.ch/work/p/pandolf/geant4/ForwardCaloUpgrade/GEANT/batchOutput_" + str(options.batchName)
if ( options.energy < 1. ) :
  dir = dir + "/" + str((int)(1000*options.energy)) + "MeV/"
else :
  dir = dir + "/" + str((int)(options.energy)) + "GeV/"


print("Creating directory: "+dir)
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/inputfiles/")
os.system("mkdir -p "+dir+"/src/")
os.system("mkdir -p "+dir+"/rootfiles/temp/")
os.system("mkdir -p "+dir+"/tables/")


pwd = os.environ['PWD']


# convert parser strings to float
act = float(options.activeLayerThickness)
abs = float(options.absorberLayerThickness)
trasv = float(options.transverseCellSize)
impact = float(options.impactPosition)

# convert floats to usable strings:
act_str = get_string(act)
abs_str = get_string(act)
trasv_str = get_string(act)
impact_str = get_string(act)



suffix = "e" + str(options.energy) + "_n"+ str(options.nLayers) + "_act" + act_str + "_abs" + abs_str + "_trasv" + trasv_str
#suffix = "e" + str(options.energy) + "_n"+ str(options.nLayers) + "_act" + str((int)(options.activeLayerThickness)) + "_abs" + str((int)(options.absorberLayerThickness)) + "_trasv" + str((int)(options.transverseCellSize))
if (int(options.impactPosition) > 0) :
    suffix = suffix + "_impact" + impact_str



# prepare the script to run
outputname = dir+"/src/submit_" + suffix + ".src"
outputfile = open(outputname,'w')
outputfile.write('#!/bin/bash\n')
outputfile.write('export LANGUAGE=C\n')
outputfile.write('export LC_ALL=C\n')
outputfile.write('cd /afs/cern.ch/work/p/pandolf/geant4/ForwardCaloUpgrade/GEANT/\n')
outputfile.write('source newconfig.sh\n')
outputfile.write("bash "+ pwd + "/scripts/simulation_batch.sh full " +str(options.energy) + " " + str(options.nLayers)+ " " + str(options.activeLayerThickness) + " " + str(options.absorberLayerThickness) + " " + str(options.transverseCellSize) + " " + str(options.impactPosition) + "\n")
outputfile.write('ls rootfiles/samplinghistos_' + suffix + '*.root | xargs -i mv {} '+dir+'/{}\n') 
outputfile.write('ls rootfiles/temp/temp_' + suffix + '*.root | xargs -i mv {} '+dir+'/{}\n') 
outputfile.write('ls tables/parameters_' + suffix + '*.dat | xargs -i mv {} '+dir+'/{}\n') 
outputfile.write('ls inputfiles/tempinput_' + suffix + '*.in | xargs -i mv {} '+dir+'/{}\n') 
outputfile.close
os.system("echo bsub -q "+options.queue+" -o "+dir+"/log/log_"+suffix+".log source "+outputname)
os.system("bsub -q "+options.queue+" -o "+dir+"/log/log_"+suffix+".log source "+outputname+" -copyInput= "+suffix)
