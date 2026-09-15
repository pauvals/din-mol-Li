#!/bin/bash

cd gcmc
rm -f Li.xyz data_gcmc

echo "Running Brownian hard-sphere simulation with GCMC reservoir"
time ../../src/dana > data_gcmc

if diff <(tail -n 596 Li.xyz) ref.xyz > /dev/null; then
  echo "Ok_GCMC"
else
  echo "FAIL_GCMC"
  exit 1
fi
