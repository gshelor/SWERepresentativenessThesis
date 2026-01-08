#!/bin/bash -l
#SBATCH --job-name=SWETreeshap		# The job name
#SBATCH --ntasks=1			# One task will be run
#SBATCH --cpus-per-task=16		# The task requires 16 CPUs
#SBATCH --hint=compute_bound		# Hyperthreading is enabled, run only on cores 
					#   - 32 total CPUs will be allocated for this job.
#SBATCH --mem-per-cpu=5000M		# Allocate 5GB of memory per CPU
					#   - 160GB = 32*5GB of memory will be allocated 
#SBATCH --time=0-08:00:00		# Max run time is 8 hours
					#   - Use --time as command line argument to override
#SBATCH --partition=cpu			# Submit job to the cpu partition
#SBATCH --mail-type=ALL			# Send mail on all state changes
#SBATCH --output=SWETreeshap.%j.out	# The output file name
#SBATCH --error=SWETreeshap.%j.err		# The error file name

#--------------------------------------------------------------------------------------------------
#
# Name: SWETreeshap
#
# Purpose: Run a containerized R job.
#
# Input: This batch script takes an R input file as an argument.
#
# Output: Output files will be written to the working directory from which the batch script was
#         submitted. The stdout file name ends in *.out. The stderr file name ends in *.err.
#         Both contain the job number.
#
# Notes: This submission script is optimized for 16 core, 160GB memory jobs. Please ensure your
#        R input file matches these dimensions. Dimensions may, optionally, be overriden
#        from the command line. See the example below.
#
# Example submission:
#
# $ sbatch SWETreeshap.sl SWE_SHAP_Slurm.R
#
# Example submission that changes maximum run time and sets a notification email address:
#
# $ sbatch --time=10-00:00:00 --mail-user=gshelor@unr.edu SWETreeshap.sl SWE_SHAP_Slurm.R
#
#--------------------------------------------------------------------------------------------------

# Load the unr-rc module to gain access to the singularity container module
module load unr-rc
# Load the singularity container module
module load singularity

# Error if the R input file does not exist or was not specified. Check stderr file for 
# error.
[[ -f ${1} ]] || { echo "R input file does not exist" >&2; exit 1; }

# Run R from the system singularity container with the R input file specified from
# the command line.
srun singularity exec /apps/R/R_tidyverse_4.4.1.sif r < ${1}