#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

## Script parameters set from environment variables
# Generate output directory (adapt this to your setup)
output=${1:-"results/test-ee"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-ee.csv

# Iterate over backends
for backend in JetReconstruction Fastjet; do
    # Durham has no R or p parameter (or implicitly R=4, p=1)
    algorithm=Durham
    cmd="julia --project src/benchmark.jl --code $backend -A $algorithm -R 4.0 -p 1.0 -m 16 --results $output $inputs"
    echo "Benchmark: $cmd"
    $cmd
    for radius in 0.2 0.4 1.0 1.5 2.0 4.0; do
        algorithm=EEKt
        for power in -1.0 0.0 1.0; do
            cmd="julia --project src/benchmark.jl --code $backend -R $radius -A $algorithm -p $power -m 16 --results $output $inputs"
            echo "Benchmark: $cmd"
            $cmd
        done
    done
done
