set key left

plot \
"band.dat" pt 6 ps 2 t "UHFk",\
-2*cos(2*pi*x/16) w l lc -1
pause -1
