#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J x.runtime.8bit
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 30:00

# Set outer loop counter.
outer=1             

# Beginning of outer loop.
for a in 64 128 256 512 -1
do
  echo "Band $outer in outer loop."
  echo "---------------------"

  # Reset inner loop counter.
  inner=1           

  # ===============================================
  # Beginning of inner loop.
  for i in 1 2 3 4 5 6 7 8 9 10
  do
    echo "Iter $inner in inner loop."
    # Run the program.  
    srun -n 1 -c 1 --cpu_bind=cores ./test-runtime-exe 100 $a
    # Increment inner loop counter.  
    let "inner+=1"  
  done
  # End of inner loop.
  # ===============================================

  let "outer+=1"    # Increment outer loop counter. 
  echo              # Space between output blocks in pass of outer loop.
done               
# End of outer loop.