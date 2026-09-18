#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
source "$repo/benchmarks/serafin/load_environment.sh"

command -v caf
command -v cafrun
caf --version

cd "$repo"
autoreconf -fi
if [[ -f Makefile ]]; then
  make distclean || true
fi
FC=caf ./configure --enable-coarray=compiler
make -j "${DANA_BUILD_JOBS:-8}"

if ! ldd src/dana | grep -q 'libcaf_mpi'; then
  echo "The executable does not appear to be linked with libcaf_mpi" >&2
  exit 1
fi

echo "Built $repo/src/dana with OpenCoarrays/MPI from the existing Spack environment"
