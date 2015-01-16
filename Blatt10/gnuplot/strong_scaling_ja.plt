set title "Strong scaling: Total runtime of partdiff with the jacobi algorithm"
set xlabel "Amount of nodes"
set ylabel "Time [t]"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
unset key
plot 'WEAK_SCALING_JA.dat' using 1:4 with points pt 7 #,  'WEAK_SCALING_GA.dat' using 1:4 with points pt 7 