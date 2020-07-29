#!/bin/bash

# Will run the following programs
### Summit_Interp.py
### Correct_Reconstruct.py 
### Summit_Climbs.py 
### getERA_Composite.py 
### Low_Extremes.py 
### Process_CMIP.py 
### Temp_Changes.py 

# From this folder:
code_dir="/home/lunet/gytm3/Everest2019/Research/OneEarth/Code/Public"

# Do it
python ${code_dir}Summit_Interp.py
python ${code_dir}Correct_Reconstruct.py
python ${code_dir}getERA_Composite.py
python ${code_dir}Low_Extremes.py
python ${code_dir}Process_CMIP.py
python ${code_dir}Temp_Changes.py
