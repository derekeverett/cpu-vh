export OMP_NUM_THREADS=$1
set OMP_NUM_THREADS $1

echo "Number Threads : ${OMP_NUM_THREADS}"
./cpu-vh --config rhic-conf/ -o output -h
