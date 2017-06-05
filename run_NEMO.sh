#!/bin/bash
#PBS -N NEMOGYRE
#PBS -l walltime=1:59:00
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -j oe

export MPI_XPMEM_ENABLED=disabled
module purge
module load intel_fc mpt hdf5 gcc
module load netcdf/4.4.1.1
module load netcdf-fortran-4.4.4_ifort


scriptdir=/home/momo/NEMOGCM/CONFIG/VENTUS_SCRIPTS
bindir=/home/momo/NEMOGCM/CONFIG/MY_GYRE_ADR/EXP00
bathydir=/home/momo/bathymetries/nemo
rm $scriptdir/NEMOGYRE.o*
cd $bindir
rm GYRE*.nc mesh*.nc *.abort*.nc
mpirun -np 1 $bindir/opa
