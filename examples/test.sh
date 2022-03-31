#!/bin/bash

# ../dana 
# diff <(tail -n 3638 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL

echo -e "Corriendo test de Ermak con pistón de partículas. Data en \e[1mdata_ermak\e[0m =)"
../dana > data_ermak 
diff <(tail -n 1206 Li.xyz) ermak_piston_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Ermak\e[0m" || echo -e "\e[1;101mFAIL_Ermak\e[0m"

echo -e "Corriendo test de Browniana con agregado de partículas. Data en \e[1mdata_brown\e[0m =)"

sed -i 's/\.true\./\.false\./g' movedor.ini

../dana > data_brown
diff <(tail -n 1206 Li.xyz) brownian_chunk_ref.xyz > /dev/null && echo -e "\e[1;103mOk_Brown\e[0m" || echo -e "\e[1;101mFAIL_Brown\e[0m"

sed -i 's/\.false\./\.true\./g' movedor.ini
