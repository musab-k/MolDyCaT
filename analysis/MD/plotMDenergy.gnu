#example gnuplot script

#required for angstroms
set encoding iso_8859_1
###Output to pdf. Will use different colours for lines:
set term pdf enhanced font "Times New Roman, 9" 
set output "diamondenergy.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
#set term postscript eps enhanced 
#set output "Example.eps"

#set nokey
#set key spacing 1.5
set autoscale
#set xrange[1:5.5]
set yrange[-8.22:-8.19]
unset log
unset label
set xtic auto
set ytic auto

#set title "All E-V curves"
set xlabel "Time (ps)"
set ylabel "Energy per atom (eV/\305)"
set style line 1 lt 1 lw 1 pt 1
set style line 2 lt 2 lw 1 pt 2
set style line 3 lt 3 lw 1 pt 3
set style line 4 lt 4 lw 1 pt 4
set style line 5 lt 5 lw 1 pt 6

f(x)=a*x+b
a=1;b=1;
fit f(x) 'diamond_energy_300K_c.dat' via a,b

plot 'diamond_energy_300K_c.dat' w l ls 1 t 'Data', f(x) ls 2 t 'Fit'
#plot 'diamond_en_300K.dat' w l ls 1 t '300 K', 'diamond_rdf_3000K.dat' w l ls 2 t '3000 K'
# removed nanotubes: 'energy_vol_nanotube.txt' w l title "nanotube", 
