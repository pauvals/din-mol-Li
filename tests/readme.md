## Para correr tests
1. Compilar en `din-mol-li` con `make` (hacer `make clean` previamente por las dudas).
1. Cambiar a esta carpeta de tests. El script `test.sh` ejecuta la regresión del
modelo conservado en `hpc-light`: dinámica Browniana de esferas rígidas,
resolución de overlaps y reservorio Gran Canónico.

En resumen, desde directorio `din-mol-li`:
		make clean
		make
		cd tests
		./test.sh
