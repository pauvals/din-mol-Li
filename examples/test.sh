
../dana 
diff <(tail -n 3638 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL
