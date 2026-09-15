## din-mol-Li
<!---
[Esp](#Esp)

[Eng](#Eng)
-->



Repo para un programa de Dinámica Molecular que permite estudiar la deposición de iones Li+ sobre un electrodo de Li metálico. 

Repo for a Molecular Dynamics program to study Li+ ions deposition onto a Li metal electrode.

### + info
There are three possible reservoir types to study, and two positions integrators.

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
causing unbounded recursion. The current implementation deliberately exchanges
a full replicated state and is a correctness baseline, not yet the final
communication pattern. Measurements show that this exchange grows poorly with
system size; compact Cartesian neighbor halos and explicit atom migration are
the next development target.

Physical atoms carry a persistent `uid` for coarray records. GEMS `gid` and
`gri` remain execution-view and group-membership indices and must not be used as
persistent distributed identities.

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
