#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 1
#SBATCH --mail-type=all
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH -t 30:00

for i in 0 1 2 3 4 5 6 7 8 9
    for x in 10 20 30 40 50 70 80 90 100
    do
        srun -n 1 -c 1 ./test $x 32 &> xavier-bmc-bw32-x$x-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for x in 10 20 30 40 50 70 80 90 100
    do
        srun -n 1 -c 1 ./test $x 64 &> xavier-bmc-bw64-x$x-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for x in 10 20 30 40 50 70 80 90 100
    do
        srun -n 1 -c 1 ./test $x 128 &> xavier-bmc-bw128-x$x-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for x in 10 20 30 40 50 70 80 90 100
    do
        srun -n 1 -c 1 ./test $x 256 &> xavier-bmc-bw256-x$x-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for x in 10 20 30 40 50 70 80 90 100
    do
        srun -n 1 -c 1 ./test $x 512 &> xavier-bmc-bw512-x$x-iter$i
    done
done