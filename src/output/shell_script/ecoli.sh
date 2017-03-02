#!/bin/bash

output=src.gpl


plotter()
{

    cat<<EOF > $output
reset
#set term postscript portrait enhanced color solid lw 2 "Helvetica" 10
#set output "|epstopdf --filter > plot.pdf"
set term pdfcairo enhanced color solid font "Arial,8" lw 2 size 7,10
set output "ecoli_aerobic_elementwise.pdf"
set xtics rotate

set multiplot
set size 0.48,0.33
set origin 0.0,0.66
set logscale y
set format x "%1.0f"
set xlabel "t [h]"
set ylabel "nr. cells / ml"
set title "E. coli"
plot "ecoli_aerobic.dat" using 1:5 w l title "implicit", "ecoli_aerobic_elementwise.dat" using 1:5 w l title "elementwise"

set origin 0.5,0.66
set format y "%1.02f"
unset logscale y
set ylabel "g/dm^{3}"
set format x "%1.0f"
set xlabel "t [h]"
set title "DOC"
plot "ecoli_aerobic.dat" using 1:4 w l notitle, "ecoli_aerobic_elementwise.dat" using 1:4 w l notitle

set size 0.48,0.32
set origin 0.0,0.33
set ylabel "mg/dm^{3}"
set xlabel "t [h]"
set format x "%1.0f"
set format y "% g"
set title "O_{2} in H_{2}O"
plot "ecoli_aerobic.dat" using 1:2 w l notitle, "ecoli_aerobic_elementwise.dat" using 1:2 w l notitle


set size 0.48,0.32
set origin 0.5,0.33
set format x "%1.0f"
set format y "% g"
set xlabel "t [h]"
set ylabel "g/dm^{3}"
set title "O_{2} in air"
plot "ecoli_aerobic.dat" using 1:3 w l notitle,  "ecoli_aerobic_elementwise.dat" using 1:3 w l notitle

set size 0.48,0.33
set origin 0.0,0.00
set format x "%1.0f"
set xlabel "t [h]"
set ylabel "nr. cells/ml"
set title "E. coli"
plot "ecoli_aerobic.dat" using 1:6 w l notitle,  "ecoli_aerobic_elementwise.dat" using 1:6 w l notitle


set size 0.48,0.33
set origin 0.5,0.00
set format x "%1.0f"
set xlabel "t [h]"
set ylabel "mg cell/ml"
set title "E. coli"
plot "ecoli_aerobic.dat" using 1:5 w l notitle, "ecoli_aerobic_elementwise.dat" using 1:5 w l notitle



unset multiplot

EOF

}


main()
{
    plotter
}

main
gnuplot $output
rm $output

exit 0;
