#!/bin/tcsh
#PBS -q class
#PBS -l nodes=1:sixcore
#PBS -l walltime=00:10:00
#PBS -N Xcelsziore
# The below chnages the working directory to the location of
# your threaded DFT program
#tO RUN INVERSE PLEASE ENSURE THAT the Transform2D 
#fucntion has a zero parameter provided as the argument
cd ThreadsTransform2D
./threadDFT2d



