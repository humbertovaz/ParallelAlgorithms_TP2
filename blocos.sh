# load do modulo


module load gnu/openmpi_eth/2.0.0
cd /home/a72095/MPI
export FILE=blocosOMP_MPI.c
mpicc $FILE -g -o blocos
mpirun -np $1 --mca btl tcp,sm,self -hostfile ${PBS_NODEFILE} ./blocos



