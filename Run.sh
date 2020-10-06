#!/bin/bash
set -e
# Will run the following programs
### Summit_Interp.py
### Correct_Reconstruct.py 
### Summit_Climbs.py 
### getERA_Composite.py 
### Low_Extremes.py 
### Process_CMIP.py 
### Temp_Changes.py 
### Changes.py

# From this folder:
code_dir="/home/lunet/gytm3/Everest2019/Research/OneEarth/Code/Public/"

# Do it
echo -e "***Starting execution...\n\n"
python ${code_dir}Summit_Interp.py
echo -e "***...Finished with Summit_Interp.py\n\n"
python ${code_dir}Correct_Reconstruct_Obs.py
echo -e "***...Finished with Correct_Reconstruct_Obs.py\n\n"
python ${code_dir}getERA_Composite.py
echo -e "***...Finished with getERA_Composite.py\n\n"
python ${code_dir}Summit_Climbs.py
echo -e "***...Finished with Summit_Climbs.py\n\n"
python ${code_dir}Low_Extremes.py
echo -e "***...Finished with Low_Extremes.py\n\n"
python ${code_dir}Process_CMIP.py
echo -e "***...Finished with Process_CMIP.py\n\n"
python ${code_dir}Temp_Changes.py
echo -e "***...Finished with Temp_Changes.py\n\n"
python ${code_dir}Changes.py
echo -e "***...Finished with Changes.py\n\n"
