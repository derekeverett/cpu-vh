export OMP_NUM_THREADS=$1
set OMP_NUM_THREADS $1

./cpu-vh --config rhic-conf/ -o output -h
