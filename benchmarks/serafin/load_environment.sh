#!/bin/bash
# Source this file from setup and Slurm jobs.

if [[ -z ${BASH_VERSION:-} ]]; then
  echo "DANA Serafin environment requires bash" >&2
  return 1 2>/dev/null || exit 1
fi

. /etc/profile
module purge

default_spack_root=/home/shared/gmm/spackrepo/spack-0.23/spack
default_spack_env=/home/shared/gmm/spackrepo/spack-0.23/serafin/v2
spack_root=${SPACK_ROOT:-$default_spack_root}
if [[ ! -f $spack_root/share/spack/setup-env.sh ]]; then
  echo "Spack was not found at $spack_root" >&2
  echo "See benchmarks/serafin/README.md for the one-time installation." >&2
  return 1 2>/dev/null || exit 1
fi

. "$spack_root/share/spack/setup-env.sh"
module load gcc spack-config

spack_env=${DANA_SPACK_ENV:-${SPACK_ENV:-$default_spack_env}}
if [[ ! -f $spack_env/spack.yaml ]]; then
  echo "Spack environment not found at $spack_env" >&2
  echo "Override it with DANA_SPACK_ENV only when needed" >&2
  return 1 2>/dev/null || exit 1
fi

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
spack_scope=$repo/benchmarks/serafin/spack-source-scope
spack -C "$spack_scope" -e "$spack_env" load opencoarrays

if ! command -v caf >/dev/null 2>&1; then
  echo "OpenCoarrays compiler wrapper 'caf' is not available" >&2
  return 1 2>/dev/null || exit 1
fi

ucx_prefix=$(spack -C "$spack_scope" -e "$spack_env" location -i ucx)
case $ucx_prefix in
  /ccad/hpcx/*)
    echo "Refusing external HPC-X UCX: $ucx_prefix" >&2
    echo "Rebuild with benchmarks/serafin/setup_open_coarrays.sh" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
