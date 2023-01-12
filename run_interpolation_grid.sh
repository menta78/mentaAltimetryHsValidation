#!/bin/bash

inFolder="/eos/jeodpp/data/projects/CLIMEX/schismRuns"
outFolder="/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/schismwwm"

for year in {1997..2015}
do
    echo $year
    rm ${outFolder}/ERA5_schismwwm*.nc
    ln -sf ${inFolder}/ERA5_schismwwm_${year}*.nc ${outFolder}/.
    envNew/bin/python main-scatter-data.py --year $year
done
