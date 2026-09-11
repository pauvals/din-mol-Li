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

Multi-image execution uses coarray intrinsics directly in `dana.F90`. The first
domain-decomposition step splits the simulation box into z-slabs.
Each image owns the `focus` atoms in its slab and keeps neighbor candidates from
its owned slab plus a z-halo of width `hs%rcut + nb_dcut`. Owned atom positions
are still synchronized through a full replicated coarray buffer after updates.
This is a correctness scaffold for neighbor/halo work, not the final scalable
communication pattern.

For multi-image runs, image 1 also builds a read-only list of unique physical
overlap pairs before the current serialized correction. Periodic ghosts are
mapped to their real atoms, pairs are assigned deterministically to the owner
of the lower global atom index, and the final profile reports local,
cross-domain, deposition, duplicate, and pair-owner counts. This diagnostic
does not change the trajectory and is the basis for distributing overlap
resolution later.

The intended production topology is a Cartesian domain decomposition. The
next topology to benchmark is an x-y process grid with domains spanning z,
which distributes the deposition surface across all images. The current
z-slab mode remains the validation baseline while ownership is generalized.

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
