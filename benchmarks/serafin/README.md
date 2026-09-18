# DANA benchmark on Serafin

Serafin provides 64 physical CPU cores and 256 GB RAM per node. These scripts
run a hybrid OpenCoarrays/MPI plus OpenMP strong-scaling benchmark on one or two
nodes. The default system tiles the 596-particle GCMC reference 24 x 24 times:
343,296 particles in a 2400 x 2400 x 100 box for 100 steps.

The larger case is intentional: the local 21,456-particle benchmark is too
small for 64 or 128 cores. Do not increase the size and image count blindly,
because the current DANA checkpoint still replicates physical atom metadata in
every image.

## 1. Transfer and unpack

From the local machine:

```sh
scp dana-hpc-96a9ae0-serafin.tar.gz \
  "$USER@serafin.ccad.unc.edu.ar:~/"
```

On Serafin:

```sh
tar -xzf ~/dana-hpc-96a9ae0-serafin.tar.gz
cd ~/dana-hpc-96a9ae0
```

## 2. Select the Spack installation and environment

The scripts accept an existing Spack installation and a dedicated environment.
For the shared Serafin Spack 0.23 tree, start from its cluster configuration
but keep DANA's roots and lock file under the user account:

```sh
export SPACK_ROOT=/home/shared/gmm/spackrepo/spack-0.23/spack
export DANA_SPACK_ENV=$HOME/spack-envs/dana

mkdir -p "$DANA_SPACK_ENV"
cp /home/shared/gmm/spackrepo/spack-0.23/serafin/v2/spack.yaml \
  "$DANA_SPACK_ENV/spack.yaml"
sed -i 's/^  specs:$/  specs: []/' "$DANA_SPACK_ENV/spack.yaml"
```

The `sed` fixes the shared manifest's empty `specs:` value, which Spack 0.23
rejects as `null`. If Spack is not installed elsewhere, the fallback remains:

```sh
git clone --branch releases/v0.21 https://github.com/spack/spack.git ~/spack
```

## 3. Build the source-only communication stack and DANA

Run, preferably inside `tmux` or `screen`:

```sh
benchmarks/serafin/setup_open_coarrays.sh
```

This adds OpenCoarrays to the selected environment, reconcretizes without
reusing installed dependencies, and builds UCX, OpenMPI, and OpenCoarrays with
GCC 12.3. The repository scope in `spack-source-scope/packages.yaml` empties
the external UCX/OpenMPI lists. It is applied by both setup and batch jobs.

This is required because Serafin's external HPC-X UCX advertises `-liberty`
and build-host paths such as `/hpc/local/oss/binutils/2.37` in `ucx-ucs.pc`.
Building OpenMPI against that metadata fails with `ld: cannot find -liberty`.
The setup aborts if concretization still selects external UCX, if DANA does not
link to `libcaf_mpi`, or if the resulting executable links under `/ccad/hpcx`.

Keep these variables exported when submitting jobs:

```sh
export SPACK_ROOT=/home/shared/gmm/spackrepo/spack-0.23/spack
export DANA_SPACK_ENV=$HOME/spack-envs/dana
```

## 4. Pilot jobs

First test one repeat, 50 steps, and a 16 x 16 system (152,576 particles):

```sh
sbatch --export=ALL,DANA_BENCH_REPEATS=1,DANA_BENCH_STEPS=50,DANA_BENCH_TILE_X=16,DANA_BENCH_TILE_Y=16 \
  benchmarks/serafin/serafin-1node.sbatch

sbatch --export=ALL,DANA_BENCH_REPEATS=1,DANA_BENCH_STEPS=50,DANA_BENCH_TILE_X=16,DANA_BENCH_TILE_Y=16 \
  benchmarks/serafin/serafin-2nodes.sbatch
```

Monitor with `squeue --me`; after a job starts, `jobmetrics JOBID` reports CPU
and memory efficiency.

## 5. Production benchmark

Submit both jobs from the repository root so they use the same executable and
input size:

```sh
sbatch benchmarks/serafin/serafin-1node.sbatch
sbatch benchmarks/serafin/serafin-2nodes.sbatch
```

One node tests `1x1`, `1x64`, `2x32`, ..., `64x1`. Two nodes test `2x64`,
`4x32`, ..., `128x1`, where each layout is `images x OpenMP threads`. Every
full-resource layout uses 64 cores per node. Results are written under
`results/serafin-JOBID-{1node,2nodes}/` with:

- `scaling.csv`: launcher and maximum per-image wall times plus residuals;
- `metadata.txt`: commit, modules, linked CAF/MPI libraries, CPU and job data;
- one raw log and timing file per repeat.

The launcher aborts if a layout exceeds 64 cores per node, cannot distribute
images evenly, or reports nonzero residual overlaps.

Combine both jobs into a median scaling table:

```sh
python3 benchmarks/serafin/summarize_scaling.py \
  results/serafin-ONE_NODE_JOB-1node/scaling.csv \
  results/serafin-TWO_NODE_JOB-2nodes/scaling.csv \
  > results/serafin-summary.csv
```

The summary reports median launcher/DANA times, maximum residuals, and speedup
relative to the `1x1` run. It rejects CSV files generated with different
particle counts or step counts.

## Size and layout overrides

Use the same tile and step values for one- and two-node strong scaling:

```sh
sbatch --export=ALL,DANA_BENCH_TILE_X=32,DANA_BENCH_TILE_Y=32,DANA_BENCH_STEPS=100 \
  benchmarks/serafin/serafin-2nodes.sbatch
```

A 32 x 32 case has 610,304 particles. Run the pilot first and inspect memory
before using it. Other useful overrides are `DANA_BENCH_REPEATS`,
`DANA_BENCH_CONFIGS`, `DANA_BENCH_GCMC_ATTEMPTS`, and `DANA_BENCH_DIR`.

The scripts use `srun --mpi=pmix`, CPU-core binding, `OMP_PLACES=cores`, and
`OMP_PROC_BIND=close`. They do not use `GFORTRAN_NUM_IMAGES`: with
OpenCoarrays/MPI, each Slurm task is one Fortran image.

## UNC references

- [Serafin hardware and software](https://wiki.ccad.unc.edu.ar/infra/clusters.html)
- [Slurm jobs and Serafin examples](https://wiki.ccad.unc.edu.ar/empezar/slurm.html)
- [Slurm partition limits](https://wiki.ccad.unc.edu.ar/infra/slurm-particiones.html)
- [Installing user software with Spack](https://wiki.ccad.unc.edu.ar/tutoriales/spack.html)
- [CPU and memory monitoring with jobmetrics](https://wiki.ccad.unc.edu.ar/tutoriales/check-efficiency.html)
