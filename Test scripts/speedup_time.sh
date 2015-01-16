#!/bin/sh
#PBS -N poisson-speedup
#PBS -o $PBS_JOBNAME.$PBS_JOBID.out
#PBS -e $PBS_JOBNAME.$PBS_JOBID.err
#PBS -q hpcintro
#PBS -l walltime=1:30:00
#PBS -l feature=XeonE5-2680
#PBS -l nodes=1:ppn=16

cd $PBS_O_WORKDIR

LOGEXT=dat

rm -f speedup_time.$LOGEXT

for thread in {1..16}
do
	echo "$thread" >> speedup_time.$LOGEXT
	time (OMP_NUM_THREADS=$thread ../poisson omp 2000 10000 0.0) 2>> speedup_time.$LOGEXT
done

exit 0
