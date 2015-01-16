set title "Weak scaling (Speedup): Total runtime of partdiff with the jacobi algorithm"
set xlabel "Amount of processes"
set ylabel "total Speedup"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
unset key
plot 'WEAK_SCALING_JA_SPEEDUP.dat' using 1:5 with points pt 7 #,  'WEAK_SCALING_GS_SPEEDUP.dat' using 1:5 with points pt 7 