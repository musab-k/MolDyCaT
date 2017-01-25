prefac=(10E3)/(2*pi)
m=1244.82
a=DUMMYA

f(x)=b*x**2+c
b=1; c=1;
fit f(x) 'pho_results_DUMMYPHO.txt' u 1:($2/108) via b,c

upp=b*2/108

wgamma=prefac*sqrt(upp/(6*m*a*a))
print wgamma 
