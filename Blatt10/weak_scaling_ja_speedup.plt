set title "Weak scaling: Speed-up of partdiff-par with the jacobi algorithm"
set xlabel "Number of processes"
set ylabel "total speed-up"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
unset key
plot 'WEAK_SCALING_JA_SPEEDUP.dat' using 1:5 with points pt 7 #,  'WEAK_SCALING_GS_SPEEDUP.dat' using 1:5 with points pt 7 
