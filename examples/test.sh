
../dana 
diff <(tail -n 2722 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL
