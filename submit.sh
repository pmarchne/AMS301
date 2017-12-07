#!/bin/bash
#
#  Name of the job (used to build the name of the standard output stream)
#$ -N N1
#
#  Number of MPI task requested
#$ -pe orte 1
#
#  The job is located in the current working directory
#$ -cwd
#
#  Merge standard error and standard output streams
#$ -j y
#
# Arguments : alpha, BC, Solver, Maxit

time mpirun -display-map ./solver 1 0 0 1000;
