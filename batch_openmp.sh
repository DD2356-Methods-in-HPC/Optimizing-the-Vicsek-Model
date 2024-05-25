
#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Second job allocation
#SBATCH -J vicsekJob_openMP
#SBATCH -t 1:00:00
#SBATCH -A edu24.DD2356
#SBATCH -p main 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -e error_file_openMP.e

# Function to run the executable multiple times
run_executable() {
    local mode=$1
    local num_iterations=$2
    local output_file=$3

    for i in $(seq 1 $num_iterations); do
        echo "Run $i - Mode: $mode" >> $output_file
        srun simulation $mode >> $output_file
        echo "" >> $output_file     # add blank line
    done
}

# Run with mode "MPI"
run_executable "openMP" 5 "output_openMP.txt"

