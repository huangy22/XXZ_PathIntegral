#set term post eps "Times-Roman" 12 color solid enh
#set term post"Times-Roman" 12
#set output "corr_k2.eps"
set xlabel "t"
set ylabel "Quan(17)"
#set title "m2a"
#set zlabel " L"
#set logscale
#set key graph 0.85, graph 0.85
#set xtics (0.56,0.57,0.58,0.59)
#set ytics (0.000,0.005,0.010,0.015)
#set xrange [0:1.0]
#set yrange [0:]
#set xrange [0.556067:0.557669]
#set nokey
set key right bottom
#set logscale y
#set logscale x
plot "./17_0.dat" u 2 with line,\
     "./17_1.dat" u 2 with line,\
     "./17_2.dat" u 2 with line,\
     "./17_3.dat" u 2 with line,\
     "./17_4.dat" u 2 with line,\
     "./17_5.dat" u 2 with line
#smooth csplines notitle  w linespoints
pause -1
