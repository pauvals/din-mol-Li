## din-mol-Li
<!---
[Esp](#Esp)

[Eng](#Eng)
-->



Repo para simular la deposición mesoscópica de iones Li+ sobre un electrodo de Li metálico.

The `hpc-light` branch contains only the model used for current HPC work:
Brownian hard-sphere dynamics, overlap resolution, and a GCMC reservoir. Ermak
dynamics, Lennard-Jones forces, piston/chunk reservoirs, constant-density
control, and their observables were intentionally removed.

### Input and output

`entrada.ini` contains the seed, deposition probability, time controls, box,
SC/SEI diffusion coefficients, neighbor skin, and initial-condition flags.
`movedor.ini` contains one line with GCMC activity and attempts per step:

    0.000627 10  ! activity, attempts

The retained outputs are `Li.xyz`, `try.dat`, and `depo.dat`.

### Compiling and testing

    # After changes in Makefile.am or configure.ac
    autoreconf -fi
    # To create Makefile
    ./configure

    make clean
    make # or make dist, distributable
    cd tests
    ./test.sh

### Coarrays

The HPC branch assumes Fortran coarray support. The default build uses GNU
Fortran's single-image coarray mode, which is the serial development path:

    ./configure
    make

This is equivalent to:

    ./configure --enable-coarray=single
    make

For a coarray runtime provided by the compiler wrapper, for example OpenCoarrays
`caf` from Spack:

    FC=caf ./configure --enable-coarray=compiler
    make

For GNU Fortran with an explicit runtime library:

    ./configure --enable-coarray=lib --with-coarray-lib=-lcaf_shmem
    make

Multi-image execution uses coarray intrinsics directly in `dana.F90`. Domain
ownership is represented by a three-dimensional Cartesian image grid. The
default topology `[1,1,num_images()]` splits the box into z slabs. Each image
owns the `focus` atoms in its domain and keeps neighboring candidates in a halo
of width `hs%rcut + nb_dcut`.

Select x-y columns at runtime with:

    DANA_DOMAIN_TOPOLOGY=xy GFORTRAN_NUM_IMAGES=4 ./dana

The image grid is factored automatically according to the lateral box aspect
ratio. Omitting `DANA_DOMAIN_TOPOLOGY` keeps z slabs.

DANA also uses OpenMP inside each image. Set `OMP_NUM_THREADS` explicitly for
multi-image benchmarks so the product of images and threads matches the CPU
allocation. For example, on a six-core allocation:

    OMP_NUM_THREADS=3 DANA_DOMAIN_TOPOLOGY=xy \
      GFORTRAN_NUM_IMAGES=2 ./dana

Multi-image overlap correction uses bounded synchronized rounds. Each image
detects conflicts involving its owned focus atoms, chooses a rollback target
from the displacement since the start of the step, applies changes only to
owned atoms, and then exchanges state before checking residual overlaps. A
conflict with no remaining displacement is reported as residual instead of
causing unbounded recursion. Routine synchronization exchanges compact boundary
halos and migration records; output checkpoints request a complete snapshot.
The current implementation still polls every image rather than only Cartesian
neighbors.

GCMC is distributed over the accessible interior of each owned reservoir
patch. A safety strip of `hs%rcut/2` is excluded at every inter-image face, so
simultaneous trial centers from neighboring patches cannot overlap. The global
`nadj` budget is apportioned by accessible volume with a largest-remainder rule;
therefore increasing the image count does not multiply the total attempt rate.
Each patch uses its accessible volume and eligible owned-particle count in the
insertion/deletion acceptance ratios.

Images first generate accepted operations against a virtual local state. They
then exchange compact ordered create/delete records, assign new UIDs by image
prefix, and replay the same record stream. The current checkpoint still keeps
the physical atom set replicated for compatibility with existing synchronization,
but it no longer broadcasts complete state for each GCMC decision. Run output
reports patch volume, accessible volume, and local attempt count. These values
must be monitored under strong scaling because fixed safety margins occupy a
larger fraction of smaller patches.

Multi-image random streams are seeded per image. Parallel trajectories are not
expected to match the serial trajectory; reservoir density, concentration
profiles, acceptance statistics, and deposition flux require ensemble-level
validation.

Physical atoms carry a persistent `uid` for coarray records. GEMS `gid(sys)` is
used only to index local owner/visibility arrays, and `gri` is used for local
group membership. Neither is a persistent distributed identity.

To run a simulation one needs `dana`, `entrada.ini` and `movedor.ini`
in the same folder.

<!---
Agregar sobre:
Contacto
Licencia
Salida de datos
dana's name has its roots in WtNV's Dana Cardinal :)
-->

Made with :metal: :mate: :battery:
