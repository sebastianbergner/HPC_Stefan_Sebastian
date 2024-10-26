#!/bin/bash

# Load all required modules
module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# compile
make --no-print-directory

REP=3
PROBLEMS="192 384 768"
RANKS="2 6 12 24 48 96"

previous_job_id=""

submit_job() {
    local ntasks=$1
    local dependency=$2
    local command=$3

    if [ -z "$dependency" ]; then
        # No dependency
        job_id=$(sbatch --ntasks=$ntasks job.sh "$command" | awk '{print $4}')
    else
        # Submit job with dependency
        job_id=$(sbatch --dependency=afterok:$dependency --ntasks=$ntasks job.sh "$command" | awk '{print $4}')
    fi

    echo $job_id
}

for i in $(seq 1 $REP); do
    for N in $PROBLEMS; do
        seq_cmd="./build/heat_stencil_2D_seq $N"
        previous_job_id=$(submit_job 1 "$previous_job_id" "$seq_cmd")

        seq_cmd="./build/heat_stencil_2D_seq_modified_access $N"
        previous_job_id=$(submit_job 1 "$previous_job_id" "$seq_cmd")
        
        for P in $RANKS; do
            par_cmd="mpiexec -n $P build/heat_stencil_2D_par_blocking $N"
            previous_job_id=$(submit_job $P "$previous_job_id" "$par_cmd")
        done
        for P in $RANKS; do
            par_cmd="mpiexec -n $P build/heat_stencil_2D_par_blockingSlices $N"
            previous_job_id=$(submit_job $P "$previous_job_id" "$par_cmd")
        done
        for P in $RANKS; do
            par_cmd="mpiexec -n $P build/heat_stencil_2D_par_blockingSlices_modified_access $N"
            previous_job_id=$(submit_job $P "$previous_job_id" "$par_cmd")
        done
    done
done
