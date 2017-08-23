LW = 2
LW2 = 1.5
LW3 = 2
LW4 = 0.5
PS = 2
EP = 50


#pointype:
#4,6,8,... is blank
#1,2,3 are lines or stars
#5,7,9,... are full

#up to 10 lline types and 12 point types
#colours, see http://www.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm

#for parameter estimation
# linetype thickens/linetype/color
# pointtype color inside /shape of point/


#line styles for solution of ODE's
set style line 1 lt 1 lw LW pi EP pt 859 ps 0.5
set style line 2 lt 2 lw LW pi EP pt 869 ps 0.5
set style line 3 lt 3 lw LW pi EP pt 879 ps 0.5
set style line 4 lt 4 lw LW pi EP pt 889 ps 0.5
set style line 5 lt 5 lw LW pi EP pt 899 ps 0.5
set style line 6 lt 6 lw LW pi EP pt 859 ps 0.5
set style line 7 lt 7 lw LW pi EP pt 869 ps 0.5
set style line 8 lt 8 lw LW pi EP pt 879 ps 0.5

#linestyle for measured data
set style line 11 lt 1 lw LW pt 859 pi EP ps PS
set style line 12 lt 2 lw LW pt 869 pi EP ps PS
set style line 13 lt 3 lw LW pt 879 pi EP ps PS
set style line 14 lt 4 lw LW pt 889 pi EP ps PS
set style line 15 lt 5 lw LW pt 899 pi EP ps PS
set style line 16 lt 6 lw LW pt 859 pi EP ps PS
set style line 17 lt 7 lw LW pt 869 pi EP ps PS
set style line 18 lt 8 lw LW pt 879 pi EP ps PS

#linestyle for results in porous media
set style line 21 lt 2 lw LW2
set style line 22 lt 4 lw LW2
set style line 23 lt 3 lw LW2
set style line 24 lt 5 lw LW2
set style line 25 lt 8 lw LW2
set style line 26 lt 12 lw LW2

#linestyle for saturation
set style line 29 lt 7 lw LW3

#linestyle for fluorescence intensity
set style line 30 lt 1 lw LW3 pt 869

#linestyle for fluorescence intensity
set style line 40 lt -1 lw LW4 pt 869


set terminal fig color big landscape metric pointsmax 10000 thickness 2 fontsize 14 size 17, 10 textspecial
set key top left spacing 1.2
set format x "%1.0f"
set xlabel "time [\\si{\\hour}]"
set samples 10000
set tics scale 0.5

solution_output= filename.".sol";
data_output= filename."_data.dat";
