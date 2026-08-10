#!/bin/bash

#parallel environment request

#$ -pe mpi_16 16

# our Job name 
#$ -N Wen.ZS

#$ -q cpl3.q
#$ -S /bin/bash
#$ -l h=(cpl2-09-ib)

#$ -cwd

cd ../

source $HOME/OpenFOAM/OpenFOAM-8.x/etc/bashrc WM_COMPILER_TYPE=ThirdParty WM_COMPILER=Gcc48 WM_LABEL_SIZE=64 WM_MPLIB=OPENMPI FOAMY_HEX_MESH=yes

export MPI_EXEC=~/OpenFOAM/ThirdParty-8.x/platforms/linux64Gcc48/openmpi-2.1.1/bin/mpirun

solver=porousPBRFoam

$MPI_EXEC -np $NSLOTS $solver -parallel > run_out.log

