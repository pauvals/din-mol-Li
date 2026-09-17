# DANA HPC benchmark

This benchmark tiles the tracked GCMC reference state into a 600 x 600 x 100
box (21,456 initial particles), runs 100 Brownian/GCMC steps, and compares
hybrid coarray/OpenMP layouts using six CPU cores.

Build with the GCC shared-memory coarray runtime, then run:

```sh
./configure --enable-coarray=lib --with-coarray-lib=-lcaf_shmem
make -j
benchmarks/hpc/run_scaling.sh
```

The GCC `libcaf_shmem` runtime used here launches its own images. Set image
count with `GFORTRAN_NUM_IMAGES`; do not wrap the executable in `oshrun`, since
that starts multiple independent coarray jobs. The script tests `1x6`, `2x3`,
`3x2`, and `6x1` (images x OpenMP threads) three times and writes a CSV plus
raw logs to a temporary directory.

Environment overrides:

```sh
DANA_BENCH_REPEATS=1 DANA_BENCH_CONFIGS="1:1 3:2" \
  DANA_BENCH_DIR=/tmp/dana-bench benchmarks/hpc/run_scaling.sh
```

Use `DANA_DOMAIN_TOPOLOGY=xy` for this wide, shallow box. The benchmark checks
the reported residual-overlap count, but it is a performance regression case,
not a substitute for the statistical physical-validation tasks in `TODO.md`.
