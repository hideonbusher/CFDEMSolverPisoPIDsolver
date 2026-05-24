#!/bin/bash

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cd $casePath
rm -f cfdem.e* cfdem.o* log*
rm -f cfdem.*
rm -f slurm* 
rm -f time*
cd $casePath/CFD
rm -rf post* cfd.mat log.liggghts core*
rm -rf [0-9].* [1-9]* VTK
rm -rf processor* clockData
rm -f *.btr
rm PID*
cd $casePath/CFD/constant/polyMesh
rm -f boundary faces neighbour owner points

cd $casePath/DEM
rm -f log.liggghts post/*
