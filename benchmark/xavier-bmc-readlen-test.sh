#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 1
#SBATCH --mail-type=all
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH -t 30:00

for i in 0 1 2 3 4 5 6 7 8 9
    for l in 500 1000 5000 10000 15000 25000 50000 100000
    do
        srun -n 1 -c 1 ./test 80 32 &> xavier-bmc-bw32-x80-len$l-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for l in 500 1000 5000 10000 15000 25000 50000 100000
    do
        srun -n 1 -c 1 ./test 80 64 &> xavier-bmc-bw64-x80-len$l-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for l in 500 1000 5000 10000 15000 25000 50000 100000
    do
        srun -n 1 -c 1 ./test 80 128 &> xavier-bmc-bw128-x80-len$l-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for l in 500 1000 5000 10000 15000 25000 50000 100000
    do
        srun -n 1 -c 1 ./test 80 256 &> xavier-bmc-bw256-x80-len$l-iter$i
    done
done

for i in 0 1 2 3 4 5 6 7 8 9
    for l in 500 1000 5000 10000 15000 25000 50000 100000
    do
        srun -n 1 -c 1 ./test 80 512 &> xavier-bmc-bw512-x80-len$l-iter$i
    done
done