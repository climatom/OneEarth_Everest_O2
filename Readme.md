# README

## Overview
This repository contains code used to assess air pressure/oxygen variability on the summit of Mount Everest. 
In the Contents section (below) the individual scripts are described. Note that raw data are not provided
here, but instructions for how to obtain the inputs are provided in Matthews et al. (submitted to One Earth). 

## Contents
### Run.sh
Calls all the python scripts, below. 
### Summit_Interp.py:
Interpolates ERA5 data to the summit of Everest (and South Col, as appropriate
### Correct_Reconstruct.py: 
Compare air pressure to the South Col, and adjust to the summit of Everest [plots Fig. 1]
### Summit_Climbs.py: 
Cross references the pressure reconstruction to the climbing history [plots Fig. 2] 
### getERA_Composite.py: 
Gets ERA5 circulation during the n most extreme low pressure events
### Low_Extremes.py: 
Takes a look at the most extreme events (lowest pressure) [plots Fig. 3]
### Process_CMIP.py: 
Calculates the sensitivity of summit pressure to global mean (near surface) air temp according to CMIP5 
### Temp_Changes.py: 
Computes observed trends in temp using HadCRUT4
### Changes.py
Computes observed trends in summit pressure and also interrogates CMIP5 results [plots Fig. 5]
