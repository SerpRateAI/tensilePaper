#!/bin/bash

# source /fp/homes01/u01/ec-benm/SerpRateAI/MicroquakesEnv/bin/activate
# source /fp/homes01/u01/ec-johnmai/.conda/envs/spec/bin

daynumber=$1
peakloc=$2

echo $daynumber
echo $peakloc

# python peakfinder.py "$daynumber" "$peakloc"
python catalogcreator.py "$daynumber"