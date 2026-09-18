#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
spack_scope=$repo/benchmarks/serafin/spack-source-scope

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

spack_env=${DANA_SPACK_ENV:-${SPACK_ENV:-}}
if [[ -z $spack_env || ! -f $spack_env/spack.yaml ]]; then
  cat >&2 <<EOF
Set DANA_SPACK_ENV to a writable, dedicated Spack environment directory.
The shared Serafin environment may be copied first; see the README.
EOF
  exit 1
fi

spack_cmd=(spack -C "$spack_scope" -e "$spack_env")
spack_spec='opencoarrays%gcc@12.3.0 ^openmpi@4.1.7%gcc@12.3.0 ^ucx~backtrace_detail'

"${spack_cmd[@]}" compiler find
"${spack_cmd[@]}" add "$spack_spec"
"${spack_cmd[@]}" concretize -f -U

spec_output=$("${spack_cmd[@]}" spec -Il opencoarrays)
printf '%s\n' "$spec_output"
if grep -Eq '\[e\].*ucx@' <<<"$spec_output"; then
  echo "Spack selected an external UCX despite the DANA source-only scope" >&2
  exit 1
fi

"${spack_cmd[@]}" install --fail-fast
"${spack_cmd[@]}" load opencoarrays

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

ucx_prefix=$("${spack_cmd[@]}" location -i ucx)
case $ucx_prefix in
  /ccad/hpcx/*)
    echo "DANA was built with the broken external HPC-X UCX: $ucx_prefix" >&2
    exit 1
    ;;
esac

if ldd src/dana | grep -q '/ccad/hpcx/'; then
  echo "DANA still links to an HPC-X library; refusing the mixed stack" >&2
  exit 1
fi

echo "Built $repo/src/dana with a source-built OpenCoarrays/OpenMPI/UCX stack"
