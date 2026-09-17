#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)

. /etc/profile
module purge

spack_root=${SPACK_ROOT:-$HOME/spack}
if [[ ! -f $spack_root/share/spack/setup-env.sh ]]; then
  cat >&2 <<EOF
Spack was not found at $spack_root.
Install it once following the UNC guide, for example:
  git clone --branch releases/v0.21 https://github.com/spack/spack.git $spack_root
EOF
  exit 1
fi

. "$spack_root/share/spack/setup-env.sh"
module load gcc spack-config

spack compiler find
spack install --fail-fast opencoarrays ^openmpi
spack load opencoarrays

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

echo "Built $repo/src/dana with OpenCoarrays/MPI"
