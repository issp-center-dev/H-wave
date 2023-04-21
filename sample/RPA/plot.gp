set dgrid3d 32,32

set term post epsf color solid enhanced

set output 'test_fig_chi0.eps'

splot \
'x' u 1:2:3 w surface lt 1 t '{/Symbol c}_0', \
1/0 notitle  

#pause -1

set output 'test_fig_chic.eps'

splot \
'x' u 1:2:4 w surface lt 2 t '{/Symbol c}_c', \
1/0 notitle  

#pause -1

set output 'test_fig_chis.eps'

splot \
'x' u 1:2:5 w surface lt 3 t '{/Symbol c}_s', \
1/0 notitle  

#pause -1

