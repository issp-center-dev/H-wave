set key left

plot \
'res.dat' u 1:(abs($3)) w lp pt 5 t'n(pi,pi)', \
'res.dat' u 1:(abs($4)) w lp pt 7 t's(pi,pi)'

pause -1
