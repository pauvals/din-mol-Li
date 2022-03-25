#!/bin/bash

# ../dana 
# diff <(tail -n 3638 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL

../dana 
diff <(tail -n 1206 Li.xyz) ermak_piston_ref.xyz > /dev/null && echo Ok_Ermak || echo FAIL_Ermak

sed -i 's/true/false/g' movedor.ini

../dana 
diff <(tail -n 1206 Li.xyz) brownian_chunk_ref.xyz > /dev/null && echo Ok_Brown || echo FAIL_Brown

sed -i 's/false/true/g' movedor.ini
