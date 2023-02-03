set key left

plot [0:16] \
"band.dat" pt 6 ps 2 t "UHFk",\
-2*cos(2*pi*x/16) w l lc -1
pause -1
