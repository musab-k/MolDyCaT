#!/bin/bash

gps='X4_omegaFinder.gnu'
grs='X4_freq.sh'
res='gruResults_X4.txt'
phonon_results='pho_X4_E_x.txt'
grunScriptOrig='grunFinderGeneral.gnu'
grunScriptNew='X4_grunScript.gnu'

cell0=3.55648

x=-0.001
dx=0.001

echo "#x cell V omega" > $res 
for (( i=0; i<=2; i++ ))
    do
        cellnew=`echo $cell0 $x | awk '{print $1*(1+$2)}'`
        V=`echo $cellnew | awk '{print $1*$1*$1}'`
        
        sed "s:DUMMYA:$cellnew:" $grs > $grs.tmp2
        sed "s:DUMMYPHO:$phonon_results:" $grs.tmp2 > $grs.tmp

        sed "s:DUMMYA:$cellnew:" $gps > $gps.tmp2
        sed "s:DUMMYPHO:$phonon_results:" $gps.tmp2 > $gps.tmp

        chmod +x $grs.tmp
        ./$grs.tmp
        gnuplot $gps.tmp 2> $gps.ans
        omega=`tail -1 $gps.ans | awk '{print $1}'`

        echo $x $cellnew $V $omega >> $res

        if [ "$x" = "0" ]; then
            omega_zero=`echo $omega | awk '{print $1}'`
        fi

        x=`echo $x $dx | awk '{print $1+$2}'`
done

sed "s:FILE:$res:" $grunScriptOrig > $grunScriptNew.tmp2
sed "s:OMEGA:$omega_zero:" $grunScriptNew.tmp2 > $grunScriptNew.tmp
gnuplot $grunScriptNew.tmp
