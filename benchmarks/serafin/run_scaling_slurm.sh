#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
binary=${DANA_BINARY:-$repo/src/dana}
repeats=${DANA_BENCH_REPEATS:-3}
configs=${DANA_BENCH_CONFIGS:-"1:64 2:32 4:16 8:8 16:4 32:2 64:1"}
nodes=${DANA_BENCH_NODES:-${SLURM_NNODES:-1}}
cores_per_node=${DANA_CORES_PER_NODE:-64}
tile_x=${DANA_BENCH_TILE_X:-24}
tile_y=${DANA_BENCH_TILE_Y:-24}
steps=${DANA_BENCH_STEPS:-100}
attempts=${DANA_BENCH_GCMC_ATTEMPTS:-10}
slurm_mpi=${DANA_SLURM_MPI:-pmix}
work=${DANA_BENCH_DIR:-${SLURM_SUBMIT_DIR:-$repo}/results/serafin-${SLURM_JOB_ID:-manual}-${nodes}node}

if [[ -z ${SLURM_JOB_ID:-} ]]; then
  echo "This launcher must run inside a Slurm allocation" >&2
  exit 1
fi
if [[ ! -x $binary ]]; then
  echo "DANA executable not found: $binary" >&2
  exit 1
fi
if [[ ${DANA_SKIP_RUNTIME_CHECK:-0} != 1 ]]; then
  "$repo/benchmarks/serafin/check_mpi_caf.sh" "$binary"
fi
for value in "$repeats" "$nodes" "$cores_per_node" "$tile_x" "$tile_y" "$steps" "$attempts"; do
  if [[ ! $value =~ ^[1-9][0-9]*$ ]]; then
    echo "Benchmark parameters must be positive integers (got '$value')" >&2
    exit 1
  fi
done
if (( nodes > SLURM_NNODES )); then
  echo "Requested $nodes benchmark nodes from a $SLURM_NNODES-node allocation" >&2
  exit 1
fi

mkdir -p "$work"

base_particles=$(wc -l < "$repo/tests/gcmc/ref.xyz")
particles=$((base_particles * tile_x * tile_y))
box_x=$((100 * tile_x))
box_y=$((100 * tile_y))

awk -v nx="$tile_x" -v ny="$tile_y" '
  { atom[NR]=$0 }
  END {
    print nx*ny*NR
    print "DANA Serafin tiled GCMC benchmark"
    for (ty=0; ty<ny; ty++) {
      for (tx=0; tx<nx; tx++) {
        for (i=1; i<=NR; i++) {
          split(atom[i], field)
          printf "%s %.12g %.12g %.12g\n", field[1], field[2]+100*tx, field[3]+100*ty, field[4]
        }
      }
    }
  }
' "$repo/tests/gcmc/ref.xyz" > "$work/pos_inic.xyz"

cat > "$work/entrada.ini" <<EOF
-104012    ! idum
1.0        ! deposition probability
5.e-5      ! timestep
$steps     ! steps
$steps     ! output interval
$box_x     ! box x
$box_y     ! box y
80         ! reservoir lower z
100        ! box upper z
250.       ! solution diffusion
250.       ! SEI diffusion
2.         ! Verlet skin
.false.    ! use pos_inic.xyz
.false.    ! no initial electrode
EOF
printf '0.000627 %s\n' "$attempts" > "$work/movedor.ini"

{
  echo "date=$(date --iso-8601=seconds)"
  echo "job_id=$SLURM_JOB_ID"
  echo "nodes=$nodes"
  echo "cores_per_node=$cores_per_node"
  echo "particles=$particles"
  echo "tile=${tile_x}x${tile_y}"
  echo "steps=$steps"
  echo "git_commit=$(git -C "$repo" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "binary=$binary"
  echo "caf=$(command -v caf || true)"
  srun --version | head -n 1
  module list 2>&1 || true
  ldd "$binary" | grep -E 'caf|mpi' || true
  lscpu
} > "$work/metadata.txt"

printf 'nodes,images,threads,total_cores,particles,steps,repeat,launcher_wall_seconds,dana_wall_seconds,residual_overlaps\n' \
  > "$work/scaling.csv"

for config in $configs; do
  if [[ $config != *:* ]]; then
    echo "Invalid config '$config'; expected images:threads" >&2
    exit 1
  fi
  images=${config%%:*}
  threads=${config##*:}
  if [[ ! $images =~ ^[1-9][0-9]*$ || ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid config '$config'; expected positive integers" >&2
    exit 1
  fi
  if (( images < nodes || images % nodes != 0 )); then
    echo "Config $config cannot be distributed evenly over $nodes nodes" >&2
    exit 1
  fi
  tasks_per_node=$((images / nodes))
  used_cores_per_node=$((tasks_per_node * threads))
  if (( used_cores_per_node > cores_per_node )); then
    echo "Config $config needs $used_cores_per_node cores/node; allocation has $cores_per_node" >&2
    exit 1
  fi
  total_cores=$((images * threads))

  for repeat in $(seq 1 "$repeats"); do
    log="$work/run-${nodes}n-${images}x${threads}-${repeat}.log"
    timing="$work/run-${nodes}n-${images}x${threads}-${repeat}.time"
    echo "Running nodes=$nodes images=$images threads=$threads repeat=$repeat"

    export OMP_NUM_THREADS="$threads"
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    export OMP_DYNAMIC=false
    export DANA_DOMAIN_TOPOLOGY=xy

    start_ns=$(date +%s%N)
    if (
      cd "$work"
      srun --mpi="$slurm_mpi" --nodes="$nodes" --ntasks="$images" \
        --ntasks-per-node="$tasks_per_node" --cpus-per-task="$threads" \
        --cpu-bind=cores --distribution=block:block --kill-on-bad-exit=1 \
        "$binary"
    ) > "$log" 2>&1; then
      run_status=0
    else
      run_status=$?
    fi
    end_ns=$(date +%s%N)

    launcher_wall=$(awk -v start="$start_ns" -v end="$end_ns" 'BEGIN {printf "%.6f", (end-start)/1.e9}')
    printf '%s\n' "$launcher_wall" > "$timing"
    if (( run_status != 0 )); then
      echo "Config $config repeat $repeat failed in srun with exit code $run_status after ${launcher_wall}s" >&2
      tail -n 80 "$log" >&2
      exit "$run_status"
    fi
    dana_wall=$(awk '/wall time:/ {if ($(NF-1)>max) max=$(NF-1)} END {printf "%.6f", max}' "$log")
    residual=$(awk -v images="$images" '
      images == 1 && /numero total de choques sin solucion:/ {
        value=$NF
        found=1
      }
      images > 1 && /residual observations:/ {
        value=$NF
        found=1
      }
      END {if(found) print value+0; else print -1}
    ' "$log")
    if (( residual < 0 )); then
      echo "Config $config repeat $repeat did not report the expected residual-overlap counter" >&2
      tail -n 80 "$log" >&2
      exit 1
    fi
    if (( residual != 0 )); then
      echo "Config $config repeat $repeat reported residual_overlaps=$residual" >&2
      tail -n 80 "$log" >&2
      exit 1
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
      "$nodes" "$images" "$threads" "$total_cores" "$particles" "$steps" \
      "$repeat" "$launcher_wall" "$dana_wall" "$residual" | tee -a "$work/scaling.csv"

    if [[ ${DANA_BENCH_KEEP_TRAJECTORY:-0} != 1 ]]; then
      rm -f "$work/Li.xyz" "$work/try.dat" "$work/depo.dat"
    fi
  done
done

echo "Results: $work/scaling.csv"
