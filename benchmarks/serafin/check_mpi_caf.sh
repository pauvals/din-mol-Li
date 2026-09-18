#!/bin/bash
set -euo pipefail

binary=${1:?Usage: check_mpi_caf.sh DANA_BINARY}

if [[ ! -x $binary ]]; then
  echo "DANA executable not found: $binary" >&2
  exit 1
fi
if ! command -v caf >/dev/null 2>&1; then
  echo "OpenCoarrays compiler wrapper 'caf' is not available" >&2
  exit 1
fi

caf_command=$(caf --show 2>&1)
case $caf_command in
  *libcaf_mpi* | *-lcaf_mpi*) ;;
  *)
    echo "The active caf wrapper does not link libcaf_mpi:" >&2
    echo "$caf_command" >&2
    exit 1
    ;;
esac

if ! binary_symbols=$(nm "$binary" 2>/dev/null); then
  echo "Could not inspect DANA symbols" >&2
  exit 1
fi
if ! grep -q '_gfortran_caf_' <<<"$binary_symbols"; then
  echo "DANA does not contain GFortran coarray runtime references" >&2
  exit 1
fi
unset binary_symbols

linked_libraries=$(ldd "$binary" 2>&1)
if grep -q 'libcaf_mpi' <<<"$linked_libraries"; then
  echo "OpenCoarrays runtime: dynamic libcaf_mpi"
elif grep -q 'libmpi' <<<"$linked_libraries"; then
  echo "OpenCoarrays runtime: static libcaf_mpi with dynamic MPI"
else
  echo "DANA does not resolve an MPI runtime:" >&2
  echo "$linked_libraries" >&2
  exit 1
fi
