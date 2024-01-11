

for iS in  2453 #469 405 #2517 1437 #3477 3541  1437  1501  #  469 405 
do


for epsT in 0.3 0.5
do

for yy in 0.1 0.3 0.5 0.7 0.9
do

outfARtoA='coop_G_ARtoA_broteT_RT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfBRtoA='coop_G_BRtoA_broteT_RT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfBRtoB='coop_G_BRtoB_broteT_RT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'

./broteT_RT_coop_G << EOF
$iS
$yy
$epsT
50
$listeps
$outfARtoA
$outfBRtoA
$outfBRtoB
EOF

done

done

done
