#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times New Roman, 9" 
set output "liquiddiamond.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

set nokey
#set key spacing 1.5
set autoscale
set xrange[1:4.5]
#set yrange[-10:0]
unset log
unset label
set xtic auto
set ytic auto

#set title "All E-V curves"
set xlabel "r (\305)"
set ylabel "g(r)"
set style line 1 lt 1 lw 1 pt 1
set style line 2 lt 2 lw 1 pt 2
set style line 3 lt 3 lw 1 pt 3
set style line 4 lt 4 lw 1 pt 4
set style line 5 lt 5 lw 1 pt 6

plot 'rdf_50000steps.dat' smooth csplines
# removed nanotubes: 'energy_vol_nanotube.txt' w l title "nanotube", 
