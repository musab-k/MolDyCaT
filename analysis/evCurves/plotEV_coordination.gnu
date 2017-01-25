#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times-Roman, 12" 
set output "coordination_energy.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set nokey
#set key top right
#set key spacing 1.5
set autoscale
set xrange[0:14]
set yrange[-9:-3]
unset log
unset label
set xtic auto
set ytic auto
set pointsize 4
set style line 4 lt -1 lw 3 pt 6

#set title "All E-V curves"
set xlabel "Coordination"
set ylabel "Energy (eV/atom)"

plot 'coordination.txt' u 1:2 w points ls 4
