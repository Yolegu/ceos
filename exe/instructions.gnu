 set terminal wxt size 800, 600
 set key off
 set xtics font 'Consolas,10'
 set ytics font 'Consolas,10'
 plot 'data.gnu' using 1:2 with lines
