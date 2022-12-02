set key right

set xlabel "Temperature T/t"
set ylabel "magnetic moment m_z"
plot \
"mag_U4.0.dat" u 1:2 pt 6 ps 2 t "U/t=4.0",\
"mag_U8.0.dat" u 1:2 pt 6 ps 2 t "U/t=8.0",\
"mag_U12.0.dat" u 1:2 pt 6 ps 2 t "U/t=12.0",\
"mag_U24.0.dat" u 1:2 pt 6 ps 2 t "U/t=24.0"
pause -1
