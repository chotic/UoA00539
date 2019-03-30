# UoA00539
Fractal Analysis of ASD

## Prerequisites

You'll need to have

 * parfor toolbox installed
 * eeglab toolbox version `14_1_2b` installed. Follow the instructions at https://sccn.ucsd.edu/eeglab/download.php to install eeglab. The scripts expect the `eeglab14_1_2b` to exist at the same level as `STEP_1_Processing_RS_v2.m`. You can create a symbolic link
```
ln -s <directory-to-eeglab14_1_2b> eeglab14_1_2b
``` 
if this is not the case.

## How to run the resting state (RS) script 

You'll need a RAW file, let's call it x.RAW, which contains EEG signal data (amplitudes vs time). Copy the RAW files to the `./` directory. 

To run the script, type
```
matlab -nodisplay -nojvm -nosplash < STEP_1_Processing_RS_v2.m
```
This will produce x.xlsx. Open the file to see values for each measure (with some Nans for some false nearest neighbours)
and print out 
```
check nansum:208.0924
```

## Adjusting the sample rate

The run time will increase significantly with smaller `downsampleRate` values. The default value is 100. Use any value >= 1 with 1 meaning
no downsampling. To set `downsampleRate` to 50, for instance, type command:
```
exit"
``` 

## Running under Slurm 

```
srun --hint=nomultithread --time=00:01:00 --cpus-per-task=10 --ntasks=1 --mem=1G \
   matlab -nodisplay -nojvm -nosplash -r "downsampleRate=100; STEP_1_Processing_RS_v2; exit"
```
where `cpus-per-task` is the number of threads and `mem` is the total memory. Use
```
squeue -u $USER
```
to check progress of job.
You can use the `sacct -j <jobid> --format=jobid,maxrss` command to infer the maximum memory used by a job from the `MaxRSS` column. 
For example:
```
sacct -j 8380550 --format=jobid,maxrss
       JobID     MaxRSS 
------------ ---------- 
8380550                 
8380550.0       790288K
```
indicating that the code took about 800MB. You want to set `mem` fairly tightly as this will let you slip into the queue faster. Allocate 4 times more wall clock time when you halve `downsampleRate`. The memory requirements don't increase significantly when halving `downsampleRate`. 

Alternatively you can edit and submit
```
sbatch STEP_1_Processing_RS_v2.sl
```
Note: you will need to adjust `downsampleRate`.
 
