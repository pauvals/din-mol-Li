
../dana 
diff <(tail -n 3942 Li.xyz) Li_ref.xyz > /dev/null && echo OK || echo FAIL
