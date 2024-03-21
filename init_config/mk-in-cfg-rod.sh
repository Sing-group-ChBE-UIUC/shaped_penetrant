#!/usr/bin/env bash
echo "Number of penetrants is: $1"
echo "Beads per penetrant: $2"
echo "length between penetrant: $3"
echo "Bead diameter:  $4"
Np=$1
AR=$2
len_pene=$3
Dp=$4

fcross=(0.125 0.333 0.571 1.000)
Ncross=(3 6 8 10)


for i in {0..5}; do
    mkdir -p f_${fcross[i]}/initcfg-Np${Np}-AR${AR}-Dp${Dp}-lp${len_pene}
    echo f_${fcross[i]}/initcfg-Np${Np}-AR${AR}-Dp${Dp}-lp${len_pene}
done

for i in {0..5}; do
    echo "Crosslink density = ${fcross[i]}"
    for j in {0..4}; do
        ./in-cfg-nonsph-cpp.out 60 30 ${Ncross[i]} 0 ${Np} ${AR} ${len_pene} ${Dp}
    done
done

