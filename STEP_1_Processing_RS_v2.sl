#!/bin/bash -e

# may need to adjust
#SBATCH --time=00:05:00
#SBATCH --job-name=step1RS
#SBATCH --output=step1RS-%A.output
#SBATCH --cpus-per-task=16          
#SBATCH --mem=18G
##SBATCH --hint=nomultithread

DOWNSAMPLE_RATE=50

#Set workers = 1 to disable parpool.
#Comment out to use max.
#WORKERS=1

# Avoid possible future version issues
module load MATLAB/2018b

# If running serial '-nojvm' can be added
matlab -nodisplay -nosplash -r "downsampleRate=$DOWNSAMPLE_RATE;setWorkers='$WORKERS';STEP_1_Processing_RS_v2; exit"

