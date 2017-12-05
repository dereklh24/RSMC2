#!/bin/bash

# Maximum number of CPUS currently allowed is 392
#SBATCH --ntasks=391
#SBATCH --constraint=infiniband
#SBATCH --cpus-per-task=1  
#SBATCH --job-name=svmodel
#SBATCH --mem-per-cpu=100MB

# Wall time (has to be in h:m:s format)
#SBATCH --time=24:00:00        
# What mail to be sent to the user
#SBATCH --mail-type=ALL               #Specify when mail should be sent to the user


cd `pathf $SLURM_SUBMIT_DIR` # go to the directory this was submitted from

#Need to load the correct module for MPI. I compiled Rmpi with gcc and openmpi, so this is what 
#I would use.
. /etc/profile.d/modules.sh

#module load mpi/openmpi/1.10.3-gcc-4.8
module load mpi/mvapich2/latest-gcc

srun --mpi=pmi2 R --no-save --silent -e 'source("02_working_example_svrv_mpi_sp500.R")'

