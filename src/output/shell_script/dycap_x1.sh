#!/bin/bash

output1=src1.gpl
output2=src2.gpl
output3=src3.gpl
output4=src4.gpl

plotter1()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_water_x1.eps"
set xtics rotate

set format x "%1.01f"
set xlabel "mg O_2 dm^{-3} water"
set ylabel "y [m]"

plot "cut_oxygen-0-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=0d", "cut_oxygen-518400-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=6d", "cut_oxygen-691200-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=8d", "cut_oxygen-777600-12" using (column(1)*1000):2 index 0 w l lw 3 title "T=9d"

EOF

cat<<EOF > $output2
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_air_x1.eps"
set xtics rotate
set format x "%g %%"
set xlabel "O_2 in air"
set ylabel "y [m]"

plot "cut_oxygen_air-0-12" using (column(1)/100*21):2 w l lw 3 title "T=0d", "cut_oxygen_air-518400-12"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=6d", "cut_oxygen_air-691200-12"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=8d", "cut_oxygen_air-777600-12"  using (column(1)/100.*21.):2 index 0 w l lw 3 title "T=9d"

EOF


cat<<EOF > $output3
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "bacterium_x1.eps"
set xtics rotate
set format y "%1.02f"
set format y "% g"
set xlabel "cells E.coli cm^{-3} water"
set ylabel "y [m]"

plot "cut_BAC-0-12" index 0 w l lw 3 title "T=0d", "cut_BAC-518400-12" index 0 w l lw 3 title "T=6d", "cut_BAC-691200-12" index 0 w l lw 3 title "T=8d", "cut_BAC-777600-12" index 0 w l lw 3 title "T=9d"
EOF


cat<<EOF > $output4
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "doc_x1.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g DOC dm^{-3} water"
set ylabel "y [m]"

plot "cut_DOC-0-12" index 0 w l lw 3 title "T=0d", "cut_DOC-518400-12" index 0 w l lw 3 title "T=6d", "cut_DOC-691200-12" index 0 w l lw 3 title "T=8d", "cut_DOC-777600-12" index 0 w l lw 3 title "T=9d"

EOF


}

plotter2()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_water_x2.eps"
set xtics rotate

set format x "%1.01f"
set xlabel "mg O_2 dm^{-3} water"
set ylabel "y [m]"

plot "cut_oxygen-0-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=0d", "cut_oxygen-518400-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=6d", "cut_oxygen-691200-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=8d", "cut_oxygen-777600-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=9d"

EOF

cat<<EOF > $output2
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_air_x2.eps"
set xtics rotate
set format x "%g %%"
set xlabel "O_2 in air"
set ylabel "y [m]"

plot "cut_oxygen_air-0-37"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=0d", "cut_oxygen_air-518400-37"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=6d", "cut_oxygen_air-691200-37"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=8d", "cut_oxygen_air-777600-37"  using (column(1)/100*21):2 index 0 w l lw 3 title "T=9d"

EOF


cat<<EOF > $output3
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "bacterium_x2.eps"
set xtics rotate
set format y "%1.02f"
set format y "% g"
set xlabel "cells E.coli cm^{-3} water"
set ylabel "y [m]"

plot "cut_BAC-0-37" index 0 w l lw 3 title "T=0d", "cut_BAC-518400-37" index 0 w l lw 3 title "T=6d", "cut_BAC-691200-37" index 0 w l lw 3 title "T=8d", "cut_BAC-777600-37" index 0 w l lw 3 title "T=9d"
EOF


cat<<EOF > $output4
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "doc_x2.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g DOC dm^{-3} water"
set ylabel "y [m]"

plot "cut_DOC-0-37" index 0 w l lw 3 title "T=0d", "cut_DOC-518400-37" index 0 w l lw 3 title "T=6d", "cut_DOC-691200-37" index 0 w l lw 3 title "T=8d", "cut_DOC-777600-37" index 0 w l lw 3 title "T=9d"

EOF


}


plotter3()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_water_x1_absolute.eps"
set xtics rotate
set format x "%1.01f"
set xlabel "mg O_2 dm^{-3} CF"
set ylabel "y [m]"

plot "cut_oxygen_absolute-0-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=0d", "cut_oxygen_absolute-518400-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=6d", "cut_oxygen_absolute-691200-12" using (column(1)*1000):2  index 0 w l lw 3 title "T=8d", "cut_oxygen_absolute-777600-12"  using (column(1)*1000):2 index 0 w l lw 3 title "T=9d"

EOF



cat<<EOF > $output3
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "bacterium_x1_absolute.eps"
set xtics rotate
set format y "%1.02f"
set format y "% g"
set xlabel "cells E.coli cm^{-3} CF"
set ylabel "y [m]"

plot "cut_BAC_absolute-0-12" index 0 w l lw 3 title "T=0d", "cut_BAC_absolute-518400-12" index 0 w l lw 3 title "T=6d", "cut_BAC_absolute-691200-12" index 0 w l lw 3 title "T=8d", "cut_BAC_absolute-777600-12" index 0 w l lw 3 title "T=9d"
EOF


cat<<EOF > $output4
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "doc_x1_absolute.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g DOC dm^{-3} CF"
set ylabel "y [m]"

plot "cut_DOC_absolute-0-12" index 0 w l lw 3 title "T=0d", "cut_DOC_absolute-518400-12" index 0 w l lw 3 title "T=6d", "cut_DOC_absolute-691200-12" index 0 w l lw 3 title "T=8d", "cut_DOC_absolute-777600-12" index 0 w l lw 3 title "T=9d"

EOF


}

plotter4()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "oxygen_water_x2_absolute.eps"
set xtics rotate
set format x "%1.01f"
set xlabel "mg O_2 dm^{-3} CF"
set ylabel "y [m]"

plot "cut_oxygen_absolute-0-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=0d", "cut_oxygen_absolute-518400-37" using (column(1)*1000):2  index 0 w l lw 3 title "T=6d", "cut_oxygen_absolute-691200-37"  using (column(1)*1000):2  index 0 w l lw 3 title "T=8d", "cut_oxygen_absolute-777600-37"  using (column(1)*1000):2 index 0 w l lw 3 title "T=9d"

EOF

cat<<EOF > $output2
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "saturation.eps"
set xtics rotate
set format x "%1.02f"
set xlabel "saturation"
set ylabel "y [m]"

plot "cut_saturation-0-20"  w l lw 3 title "water saturation"

EOF


cat<<EOF > $output3
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "bacterium_x2_absolute.eps"
set xtics rotate
set format y "%1.02f"
set format y "% g"
set xlabel "cells E.coli cm^{-3} CF"
set ylabel "y [m]"

plot "cut_BAC_absolute-0-37" index 0 w l lw 3 title "T=0d", "cut_BAC_absolute-518400-37" index 0 w l lw 3 title "T=6d", "cut_BAC_absolute-691200-37" index 0 w l lw 3 title "T=8d", "cut_BAC_absolute-777600-37" index 0 w l lw 3 title "T=9d"
EOF


cat<<EOF > $output4
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "doc_x2_absolute.eps"
set xtics rotate

set format x "%1.02f"
set xlabel "g DOC dm^{-3} CF"
set ylabel "y [m]"

plot "cut_DOC_absolute-0-37" index 0 w l lw 3 title "T=0d", "cut_DOC_absolute-518400-37" index 0 w l lw 3 title "T=6d", "cut_DOC_absolute-691200-37" index 0 w l lw 3 title "T=8d", "cut_DOC_absolute-777600-37" index 0 w l lw 3 title "T=9d"

EOF


}

plotter5()
{

cat<<EOF > $output1
reset
set term postscript enhanced portrait color solid lw 2 "Helvetica" 18
set output "air_diffusion.eps"
set xtics rotate
set logscale x
set format x "%g "
set xlabel "difussion"
set ylabel "y[m]"

plot "cut_D-0-20"  w l lw 3 title "air diffusion"

EOF





}


main()
    {

	plotter1
	gnuplot $output1
	gnuplot $output2
	gnuplot $output3
	gnuplot $output4
	plotter2
	gnuplot $output1
	gnuplot $output2
	gnuplot $output3
	gnuplot $output4
	plotter3
	gnuplot $output1
	gnuplot $output3
	gnuplot $output4
	plotter4
	gnuplot $output1
	gnuplot $output2
	gnuplot $output3
	gnuplot $output4
	plotter5
	gnuplot $output1

    }
    mv cut_* ./solution_cut/
    cd solution_cut

    main

    rm $output1 $output2 $output3 $output4
    mv *.eps ./../pict/
    cd ..

    exit 0;
