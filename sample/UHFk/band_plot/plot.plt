set key left bottom

plot [:] \
"band.dat" u 2:4 pt 6 ps 2 t "UHFk",\
-2*cos(pi*x) w l lc -1 t "Exact"
pause -1
