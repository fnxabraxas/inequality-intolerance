

for iS in  2453 #469 405 #2517 1437 #3477 3541  1437  1501  #  469 405 
do


for epsT in 0.3 0.5
do

for yy in 0.1 0.3 0.5 0.7 0.9
do

outfARtoA='coop_G_ARtoA_broteR_TT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfARtoB='coop_G_ARtoB_broteR_TT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfBRtoA='coop_G_BRtoA_broteR_TT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfBRtoB='coop_G_BRtoB_broteR_TT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'
outfAMtoA='coop_G_AMtoA_broteR_TT_epsT_'$epsT'_y_'$yy'_stg_'$iS'.dat'

./broteR_TT_coop_G << EOF
$iS
$yy
$epsT
50
$outfARtoA
$outfARtoB
$outfBRtoA
$outfBRtoB
$outfAMtoA
EOF

done

done

done
