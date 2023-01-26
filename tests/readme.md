## Para correr tests
1. Compilar en `din-mol-li` con `make` (hacer `make clean` previamente por las dudas).
1. Cambiar a esta carpeta de tests. El script `test.sh` corre simulaciones que estarán en carpetas
distintas; una usando algoritmo de Ermak con un pistón que inyecta partículas; otra con algoritmo 
Browniano y agregado de bloques de partículas; otra con algoritmo Browniano y reservorio Gran Canónico.

En resumen, desde directorio `din-mol-li`:
		make clean
		make
		cd tests
		./test.sh
