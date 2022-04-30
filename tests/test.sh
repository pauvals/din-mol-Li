#!/bin/bash

tests=${1:-ermak brown gcmc}
 
if [[ "$tests" == *"ermak"* ]]; then
	echo -e "Corriendo test de Ermak con pistón de partículas. Data en \e[1mermak/data_ermak\e[0m =)"
	cd ermak
	rm -f Li.xyz data_ermak
	time ( ../../dana > data_ermak )
	diff <(tail -n 1206 Li.xyz) ermak_piston_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Ermak\e[0m" || echo -e "\e[1;101mFAIL_Ermak\e[0m"
	cd ..
fi


if [[ "$tests" == *"brown"* ]]; then
	echo -e "Corriendo test de Browniana con agregado de partículas. Data en \e[1mbrown/data_brown\e[0m =)"
	cd brown
	rm -f Li.xyz data_brown

	time ( ../../dana > data_brown )
	diff <(tail -n 3012 Li.xyz) brownian_chunk_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Brown\e[0m" || echo -e "\e[1;101mFAIL_Brown\e[0m"
	cd ..
fi
 

if [[ "$tests" == *"gcmc"* ]]; then
	echo -e "Corriendo test de $1. Data en \e[1mbrown/data_$1\e[0m =)"
	cd gcmc
	rm -f Li.xyz data_$1

	time ( ../../dana > data_$1 )
	diff <(tail -n 3012 Li.xyz) ref.xyz > /dev/null && echo -e	"\e[1;103mOk_$1\e[0m" || echo -e "\e[1;101mFAIL_$1\e[0m"
	cd ..
fi
 
#sed -i 's/\.false\./\.true\./g' movedor.ini
echo "Saliendo de carpeta de tests..."
cd ..

echo "¡Listo!"
