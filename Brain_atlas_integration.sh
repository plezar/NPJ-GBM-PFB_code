#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe mpi-24 120     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N integration       # Specify job name

module load R/4.1.1/gcc

Rscript Brain_atlas_integration.R
