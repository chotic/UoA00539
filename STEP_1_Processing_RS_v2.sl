#!/bin/bash -e
#
# may need to adjust
#SBATCH --time=00:05:00

#SBATCH --job-name=step1RS
##SBATCH --error=step1RS-%A.error
#SBATCH --output=step1RS-%A.output
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2          #
#SBATCH --mem=10G
##SBATCH --hint=nomultithread

DOWNSAMPLE_RATE=50
WORKERS=1

# change to "module load MATLAB" on mahuika
module load MATLAB
matlab -nodisplay -nosplash -r "downsampleRate=$DOWNSAMPLE_RATE;manualWorkerThreads='$WORKERS';STEP_1_Processing_RS_v2; exit"

