

for iS in  2517 1437 #3477 3541  1437  1501  #  469 405 
do

listeps='list_eps_'$iS'.dat'

for yy in 0.1 0.3 0.5 0.7 0.9
do

outfARtoA='coop_ARtoA_broteT_RT_'$yy'_'$iS'.dat'
outfBRtoA='coop_BRtoA_broteT_RT_'$yy'_'$iS'.dat'
outfBRtoB='coop_BRtoB_broteT_RT_'$yy'_'$iS'.dat'

./broteT_RT_coop << EOF
$iS
$yy
0.01d0
0.99d0
50
$listeps
$outfARtoA
$outfBRtoA
$outfBRtoB
EOF

done

done
