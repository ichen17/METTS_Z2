#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=5-00:00:00   # walltime limit (HH:MM:SS) 130 for 20
#SBATCH --nodes=1   # number of nodes - I always use just 1
#SBATCH --ntasks-per-node=8    #  processor core(s) per node - in our case this is the discretization of q-points
#SBATCH --mail-user=ichen@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=speedy,swift,bigram,gpu
#SBATCH --mem=100G
#SBATCH --array=34,35
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#submit with sbatch script.py


#record start date and host
date
hostname
# unload somthing I don't need
module purge
# load python 3.8
module load julia/1.8.3-py310-tgf4scr


# run program
julia METTS_E_m_N_var_mu.jl ${SLURM_ARRAY_TASK_ID}
#os.system('python -m script_a.py')

#sbatch script.py
#squeue --user=ichen
#record end date
date
