
for b in 1.5 2 # 3
do

for iS in    2517 1437 #3477 3541  1437  1501  #  469 405 
do

outf='broteR_TT_'$b'_'$iS'_101.dat'

./broteR_TT << EOF
$iS
$b
0.01d0
0.99d0
101
0.001d0
0.999d0
0.01d0
$outf
EOF

done

done
