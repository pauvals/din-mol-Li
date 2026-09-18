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

"$repo/benchmarks/serafin/check_mpi_caf.sh" src/dana

echo "Built $repo/src/dana with OpenCoarrays/MPI from the existing Spack environment"
