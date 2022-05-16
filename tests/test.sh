#!/bin/bash

tests=${1:-ermak brown gcmc}
 
if [[ "$tests" == *"ermak"* ]]; then
	echo -e "Corriendo test de Ermak con pistón de partículas. Data en \e[1mermak/data_ermak\e[0m =)"
	cd ermak
	rm -f Li.xyz data_ermak
	time ( ../../src/dana > data_ermak )
	diff <(tail -n 1206 Li.xyz) ermak_piston_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Ermak\e[0m" || echo -e "\e[1;101mFAIL_Ermak\e[0m"
	cd ..
fi


if [[ "$tests" == *"brown"* ]]; then
	echo -e "Corriendo test de Browniana con agregado de partículas. Data en \e[1mbrown/data_brown\e[0m =)"
	cd brown
	rm -f Li.xyz data_brown

	time ( ../../src/dana > data_brown )
	diff <(tail -n 3012 Li.xyz) brownian_chunk_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Brown\e[0m" || echo -e "\e[1;101mFAIL_Brown\e[0m"
	cd ..
fi
 

if [[ "$tests" == *"gcmc"* ]]; then
  x=gcmc
	echo -e "Corriendo test de $x. Data en \e[1mgcmc/data_$x\e[0m =)"
	cd gcmc
	rm -f Li.xyz data_$x

	time ( ../../src/dana > data_$x )
	diff <(tail -n 1206 Li.xyz) ref.xyz > /dev/null && echo -e	"\e[1;103mOk_$x\e[0m" || echo -e "\e[1;101mFAIL_$x\e[0m"
	cd ..
fi
 
echo "Saliendo de carpeta del test..."
cd ..

echo "¡Listo!"
