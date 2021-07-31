
../dana 
diff <(tail -n 2002 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL
