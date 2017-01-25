#example gnuplot script

###Output to pdf. Will use different colours for lines:
#set term pdf enhanced 
#set output "Example.pdf"

###Output to eps. These can easily be included in a TeX document. They will be black and white. 
set term postscript eps enhanced 
set output "Example.eps"

#set key cent right
#set key spacing 1.5
set autoscale
set xrange[0:50]
#set yrange[0:1]
unset log
unset label
set xtic auto
set ytic auto

set title "Super cool graph"
set xlabel "Temperature ({/Symbol \260} C)"
set ylabel "Number of bananas (m^{-3})"

plot x**2 t 'Data', (1.1*x**2) t 'Expected figures' 
