f(x)=a*x+b*x**2
fit [0:0.6] f(x) "../corr_k2.txt" using ($1)**(-2.0):2 via "fit.par"
set xrange [0:1.0]
set yrange [0:]
plot "./corr_k2.txt" using ($1)**(-2.0):2, f(x)
pause -1
