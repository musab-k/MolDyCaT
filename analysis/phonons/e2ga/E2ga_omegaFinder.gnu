prefac=(1E3)/(2*pi)
m=1244.82
a=DUMMYA

f(x)=b*x**2+c
b=1; c=1;
fit f(x) 'DUMMYPHO' u 1:($2/300) via b,c

upp=b*2

wgamma=prefac*sqrt(upp/(2*m*a*a))
print wgamma 
