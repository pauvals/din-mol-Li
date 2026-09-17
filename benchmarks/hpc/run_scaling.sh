#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
binary=${DANA_BINARY:-$repo/src/dana}
repeats=${DANA_BENCH_REPEATS:-3}
configs=${DANA_BENCH_CONFIGS:-"1:6 2:3 3:2 6:1"}
work=${DANA_BENCH_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/dana-hpc.XXXXXX")}

mkdir -p "$work"

awk '
  { atom[NR]=$0 }
  END {
    print 36*NR
    print "DANA HPC tiled GCMC benchmark"
    for (ty=0; ty<6; ty++) {
      for (tx=0; tx<6; tx++) {
        for (i=1; i<=NR; i++) {
          split(atom[i], field)
          printf "%s %.12g %.12g %.12g\n", field[1], field[2]+100*tx, field[3]+100*ty, field[4]
        }
      }
    }
  }
' "$repo/tests/gcmc/ref.xyz" > "$work/pos_inic.xyz"

cat > "$work/entrada.ini" <<'EOF'
-104012    ! idum
1.0        ! deposition probability
5.e-5      ! timestep
100        ! steps
100        ! output interval
600        ! box x
600        ! box y
80         ! reservoir lower z
100        ! box upper z
250.       ! solution diffusion
250.       ! SEI diffusion
2.         ! Verlet skin
.false.    ! use pos_inic.xyz
.false.    ! no initial electrode
EOF
printf '0.000627 10\n' > "$work/movedor.ini"

printf 'images,threads,repeat,wall_seconds,residual_overlaps\n' > "$work/scaling.csv"
for config in $configs; do
  images=${config%%:*}
  threads=${config##*:}
  for repeat in $(seq 1 "$repeats"); do
    log="$work/run-${images}x${threads}-${repeat}.log"
    (
      cd "$work"
      env OMP_NUM_THREADS="$threads" \
        GFORTRAN_NUM_IMAGES="$images" \
        DANA_DOMAIN_TOPOLOGY=xy \
        "$binary"
    ) > "$log" 2>&1

    wall=$(awk '/wall time:/ {if ($(NF-1)>max) max=$(NF-1)} END {printf "%.6f", max}' "$log")
    residual=$(awk '/residual observations:/ {value=$NF} END {print value+0}' "$log")
    printf '%s,%s,%s,%s,%s\n' "$images" "$threads" "$repeat" "$wall" "$residual" \
      | tee -a "$work/scaling.csv"
  done
done

printf 'Results: %s/scaling.csv\n' "$work"
