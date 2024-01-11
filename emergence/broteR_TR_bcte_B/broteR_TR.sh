
for b in 1.5 # 2 3
do

for iS in  469 # 1437  2517 # 2453 3477 3541 1437 1501 #  469 405 
do

outf='broteR_TR_'$b'_'$iS'.dat'

./broteR_TR << EOF
$iS
$b
0.01d0
0.99d0
21
0.001d0
0.999d0
0.01d0
$outf
EOF

done

done
