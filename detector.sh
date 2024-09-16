#!/bin
# source /fp/homes01/u01/ec-benm/SerpRateAI/MicroquakesEnv/bin/activate
# source /fp/homes01/u01/ec-johnmai/.peakconda/envs/spec/bin

daynumber=$1
peakloc=$2
peakdistance=$3

echo $daynumber
echo $peakloc
echo $peakdistance

python peakfinder.py "$daynumber" "$peakloc" "$peakdistance"
# python catalogcreator.py "$daynumber"
python everything_pipeline.py "$daynumber" "$peakloc"