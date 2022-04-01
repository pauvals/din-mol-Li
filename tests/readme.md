## Para correr tests
### _En construcción_
1. Compilar en `din-mol-li` con `make` (hacer `make clean` previamente por las dudas).
1. Cambiar a esta carpeta de tests. El script `test.sh` corre dos simulaciones que estarán en carpetas
distintas; una usando algoritmo de Ermak con un pistón que inyecta partículas, y otra con algoritmo 
Browniano y agregado de bloques de partículas.

En resumen, desde directorio `din-mol-li`:
		make clean
		make
		cd tests
		./test.sh
