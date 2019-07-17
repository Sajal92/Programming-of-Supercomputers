#!/bin/bash
#@ wall_clock_limit = 00:30:00
#@ job_name = pos-gauss-mpi-intel
#@ job_type = MPICH
#@ output = out_gauss_64_intel_$(jobid).out
#@ error = out_gauss_64_intel_$(jobid).out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = gauss
#@ minimize_time_to_solution = yes
#@ notification = never
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.intel

#Uncomment the following lines for tracing
#module load scorep
#ulimit -c unlimited
#export SCOREP_TIMER="clock_gettime"
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_ENABLE_TRACING=true
#export SCOREP_TOTAL_MEMORY=2GB

date
module list

echo "Tests starting..."
echo

mpiexec -n 8 ./gauss ./ge_data/size64x64
date
mpiexec -n 16 ./gauss ./ge_data/size64x64
date
mpiexec -n 32 ./gauss ./ge_data/size64x64
date
mpiexec -n 64 ./gauss ./ge_data/size64x64
date

mpiexec -n 8 ./gauss ./ge_data/size512x512
date
mpiexec -n 16 ./gauss ./ge_data/size512x512
date
mpiexec -n 32 ./gauss ./ge_data/size512x512
date
mpiexec -n 64 ./gauss ./ge_data/size512x512
date

mpiexec -n 8 ./gauss ./ge_data/size1024x1024
date
mpiexec -n 16 ./gauss ./ge_data/size1024x1024
date
mpiexec -n 32 ./gauss ./ge_data/size1024x1024
date
mpiexec -n 64 ./gauss ./ge_data/size1024x1024
date

mpiexec -n 8 ./gauss ./ge_data/size2048x2048
date
mpiexec -n 16 ./gauss ./ge_data/size2048x2048
date
mpiexec -n 32 ./gauss ./ge_data/size2048x2048
date
mpiexec -n 64 ./gauss ./ge_data/size2048x2048
date

mpiexec -n 8 ./gauss ./ge_data/size4096x4096
date
mpiexec -n 16 ./gauss ./ge_data/size4096x4096
date
mpiexec -n 32 ./gauss ./ge_data/size4096x4096
date
mpiexec -n 64 ./gauss ./ge_data/size4096x4096
date

mpiexec -n 8 ./gauss ./ge_data/size8192x8192
date
mpiexec -n 16 ./gauss ./ge_data/size8192x8192
date
mpiexec -n 32 ./gauss ./ge_data/size8192x8192
date
mpiexec -n 64 ./gauss ./ge_data/size8192x8192
date

echo
echo "Tests ended."

