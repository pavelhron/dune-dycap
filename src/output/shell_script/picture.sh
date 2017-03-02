#!/bin/bash

output1=src1.gpl
output2=src2.gpl
output3=src3.gpl
output4=src4.gpl

plotter()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_water.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "mg O_2 in water dm^{-3} fluid"
set ylabel "y [m]"

plot "cut_oxygen-0-20" index 0 w l lw 3 title "T=0d", "cut_oxygen-518400-20" index 0 w l lw 3 title "T=6d", "cut_oxygen-172800-20" index 0 w l lw 3 title "T=8d"

EOF

cat<<EOF > $output2
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_air.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g O_2 in air dm^{-3} fluid"
set ylabel "y [m]"

plot "cut_oxygen_air-0-20" index 0 w l lw 3 title "T=0d", "cut_oxygen_air-518400-20" index 0 w l lw 3 title "T=6d", "cut_oxygen-172800-20" index 0 w l lw 3 title "T=8d"

EOF


cat<<EOF > $output3
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "bacterium.eps"
set xtics rotate
set format y "%1.02f"
set format y "% g"
set xlabel "cells E.coli dm^{-3} fluid"
set ylabel "y [m]"

plot "cut_BAC-0-20" index 0 w l lw 3 title "T=0d", "cut_BAC-518400-20" index 0 w l lw 3 title "T=6d", "cut_BAC-172800-20" index 0 w l lw 3 title "T=8d"
EOF


cat<<EOF > $output4
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "doc.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g DOC dm^{-3} fluid"
set ylabel "y [m]"

plot "cut_DOC-0-20" index 0 w l lw 3 title "T=0d", "cut_DOC-518400-20" index 0 w l lw 3 title "T=6d", "cut_DOC-172800-20" index 0 w l lw 3 title "T=8d"

EOF


}


main()
    {

	plotter
    }

    main
    gnuplot $output1
gnuplot $output2
  gnuplot $output3
gnuplot $output4
    rm $output1 $output2

cd ./pict/
rm *.eps
cd ..
mv *.eps ./pict/
cd pict/
cp *.eps ../../../../../../svn/pavel/dycap/hgs_colloquium_poster2011/pics
cd ..
mv *vtu ./vtk/
mv *pvd ./vtk/

    exit 0;
