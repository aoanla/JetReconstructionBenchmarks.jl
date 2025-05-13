#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

## Script parameters set from environment variables
# Generate output directory (adapt this to your setup)
output=${1:-"results/test-antikt-ipmi"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-subset.csv

# Iterate over backends, including Python
#for backend in JetReconstruction Fastjet AkTPython AkTNumPy; do
for backend in AkTNumPy; do
    for strategy in N2Plain N2Tiled; do
	 radius=4.0
#        for radius in 0.5 1.0 2.0 4.0; do
            # Python is slow, so do less trials
            if [ $backend == "JetReconstruction" -o $backend == "Fastjet" ]; then
                trials=60
            else
                trials=4
            fi
            algorithm=AntiKt
            cmd="julia -t 3 --project=. src/benchmark.jl --code $backend --radius $radius --algorithm $algorithm --strategy $strategy -m $trials --results $output --ipmi $inputs"
            echo "Benchmark: $cmd"
            $cmd
#        done
    done
done
