#!/bin/bash
# Source this file from setup and Slurm jobs.

if [[ -z ${BASH_VERSION:-} ]]; then
  echo "DANA Serafin environment requires bash" >&2
  return 1 2>/dev/null || exit 1
fi

. /etc/profile
module purge

spack_root=${SPACK_ROOT:-/home/shared/gmm/spackrepo/spack-0.23/spack}
spack_env=${DANA_SPACK_ENV:-/home/shared/gmm/spackrepo/spack-0.23/serafin/v2}
if [[ ! -f $spack_root/share/spack/setup-env.sh ]]; then
  echo "Spack was not found at $spack_root" >&2
  echo "See benchmarks/serafin/README.md for the one-time installation." >&2
  return 1 2>/dev/null || exit 1
fi
if [[ ! -f $spack_env/spack.yaml ]]; then
  echo "Spack environment was not found at $spack_env" >&2
  return 1 2>/dev/null || exit 1
fi

. "$spack_root/share/spack/setup-env.sh"
module load gcc spack-config
spack -e "$spack_env" load opencoarrays

if ! command -v caf >/dev/null 2>&1; then
  echo "OpenCoarrays compiler wrapper 'caf' is not available" >&2
  return 1 2>/dev/null || exit 1
fi
