#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Second job allocation
#SBATCH -J vicsekJob_openMP
#SBATCH -t 1:00:00
#SBATCH -A edu24.DD2356
#SBATCH -p main 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --nodes=1
#SBATCH -e error_file_openMP.e

# Function to run the executable multiple times
run_executable() {
    local mode=$1
    local num_iterations=$2
    local output_file=$3

    # remove old values
    > $output_file

    for i in $(seq 1 $num_iterations); do
        echo "Run $i - Mode: $mode" >> $output_file
        srun bin/bench.out $mode >> $output_file
        echo "" >> $output_file     # add blank line
    done
}

# Outer loop over number of threads = powers of 2
for num_threads in 1 2 4 8 16 32 64 128; do
    echo "Running with $num_threads threads:"
    export OMP_NUM_THREADS=$num_threads
    run_executable "openMP" 5 "output_openMP_${num_threads}.txt"
done
