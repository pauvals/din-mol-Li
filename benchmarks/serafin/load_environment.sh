#!/bin/bash
# Source this file from setup and Slurm jobs.

if [[ -z ${BASH_VERSION:-} ]]; then
  echo "DANA Serafin environment requires bash" >&2
  return 1 2>/dev/null || exit 1
fi

case $- in
  *u*) restore_nounset=1; set +u ;;
  *) restore_nounset=0 ;;
esac
. /etc/profile
if ((restore_nounset)); then
  set -u
fi
unset restore_nounset
module purge

spack_root=${SPACK_ROOT:-/home/shared/gmm/spackrepo/spack-0.23/spack}
spack_env=${DANA_SPACK_ENV:-/home/shared/gmm/spackrepo/spack-0.23/serafin/v2}
if [[ ! -x $spack_root/bin/spack ]]; then
  echo "Spack was not found at $spack_root" >&2
  echo "See benchmarks/serafin/README.md for the one-time installation." >&2
  return 1 2>/dev/null || exit 1
fi
if [[ ! -f $spack_env/spack.yaml ]]; then
  echo "Spack environment was not found at $spack_env" >&2
  return 1 2>/dev/null || exit 1
fi

module load gcc spack-config
if ! spack_load_commands=$("$spack_root/bin/spack" -e "$spack_env" load --sh opencoarrays); then
  echo "Could not load OpenCoarrays from $spack_env" >&2
  return 1 2>/dev/null || exit 1
fi
eval "$spack_load_commands"
unset spack_load_commands

if ! command -v caf >/dev/null 2>&1; then
  echo "OpenCoarrays compiler wrapper 'caf' is not available" >&2
  return 1 2>/dev/null || exit 1
fi
