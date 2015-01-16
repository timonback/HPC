set title "Communication: Total runtime of partdiff with the jacobi algorithm over multiple nodes"
set xlabel "Amount of nodes"
set ylabel "Time [t]"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
unset key
plot 'COMMUNICATION_A_JA.dat' using 2:4 with points pt 7 #,  'COMMUNICATION_A_GS.dat' using 2:4 with points pt 7 