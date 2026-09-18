# HPC TODO

## Active

- [x] Establish a reproducible strong-scaling benchmark for the reduced
  Brownian hard-sphere plus GCMC model.
- [x] Profile serial, coarray, and hybrid runs with equal CPU allocations.
- [x] Reduce redundant synchronization, global apply scans, and neighbor
  rebuilds in the overlap solver.
- [x] Record repeated wall times, speedup, residual overlaps, and overlap
  rounds for the case-study poster.
- [ ] Run and archive the one- and two-node OpenCoarrays/MPI scaling campaigns
  on Serafin, including `scaling.csv`, metadata, raw logs, and `jobmetrics`.
- [ ] Instrument bytes/records exchanged and per-step halo and migration
  counts.

## Physical validation

These checks are required before claiming that the parallel reservoir preserves
the production physics, but they do not block the current HPC implementation
and scaling work.

- [ ] Compare reservoir mean density across one, two, and four images over
  multiple independent seeds.
- [ ] Compare x, y, and z concentration profiles, including bins around
  inter-image safety margins.
- [ ] Compare GCMC insertion/deletion proposal and acceptance rates.
- [ ] Compare deposition flux and accumulated deposition with uncertainties.
- [ ] Estimate Brownian mixing time in the reservoir and compare it with the
  depletion/deposition time scale.
- [ ] Define the minimum acceptable accessible-volume fraction `Va/Vp` for
  production decompositions.

## Later architecture

- [ ] Replace all-image polling with direct Cartesian-neighbor communication.
- [ ] Remove the remaining replicated physical atom set once halo creation,
  deletion, migration, and output gathering are independently validated.
- [ ] Replace stateful per-image random streams with keyed/counter random
  values if decomposition-independent reproducibility becomes a requirement.
