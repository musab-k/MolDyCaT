#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 9" 
set output "evCurves.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set nokey
#set key spacing 1.5
set autoscale
set xrange[1.2:2.2]
set yrange[-10:0]
unset log
unset label
set xtic auto
set ytic auto

#set title "All E-V curves"
set xlabel "Nearest neighbour separation (\305)"
set ylabel "Energy (eV/atom)"
set label "BCC" at 1.83,-2 rotate by -63
set label "FCC" at 1.9, -0.5
set label "HCP" at 1.81,-1.2
set label "{/Symbol b}-Sn" at 1.62,-0.5
set label "SC" at 1.53,-1.5
set label "Diamond" at 1.65,-8.5
set label "Graphite" at 1.32,-8.7
set label "Linear Chain" at 1.41,-6.6
set style line 1 lt 1 lw 1 pt 1
set style line 2 lt 2 lw 1 pt 2
set style line 3 lt 3 lw 1 pt 3
set style line 4 lt 4 lw 1 pt 4
set style line 5 lt 5 lw 1 pt 6

plot 'energy_vol_BCC.txt' w l ls 1 title "BCC", 'energy_vol_FCC.txt' w l ls 2 title "FCC", 'energy_vol_HCP.txt' w l ls 3 title "HCP", 'energy_vol_Graphene.txt' w l title "graphite", 'energy_vol_Sn.txt' w l title "{/Symbol b}-Sn", 'energy_vol_SC.txt' w l title "SC", 'energy_vol_diamond512.txt' w l title "diamond",'energy_vol_linear_chain.txt' w l title "linear chain"
 # removed nanotubes: 'energy_vol_nanotube.txt' w l title "nanotube", 
