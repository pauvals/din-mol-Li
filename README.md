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
