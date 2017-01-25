f(x)=a*x+b
a=1; b=1;
fit f(x) 'FILE' u 3:4 via a,b
V=4.257*4.257*3.35/sqrt(3)
w=OMEGA
g=-a*V/w
print g

plot 'FILE' u 3:4, f(x)
