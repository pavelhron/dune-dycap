set pointsize 3
PT = 9
LW = 5

#pointype:
#4,6,8,... is blank
#1,2,3 are lines or stars
#5,7,9,... are full

#up to 10 lline types and 12 point types
#colours, see http://www.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm

set style line 1 lt 1 lw LW pt 4 linecolor rgb "red"
set style line 2 lt 2 lw LW pt 6 linecolor rgb "orange"
set style line 3 lt 3 lw LW pt 8 linecolor rgb "violet"
set style line 4 lt 1 lw LW pt 10 linecolor rgb "green"
set style line 5 lt 5 lw LW pt 12 linecolor rgb "blue"
set style line 6 lt 5 lw LW pt 4 linecolor rgb "black"
set style line 7 lt 5 lw LW pt 4 linecolor rgb "violet"
set style line 8 lt 5 lw LW pt 4 linecolor rgb "cyan"
set style line 9 lt 5 lw LW pt 4 linecolor rgb "#006400"

set style line 10 lt 1 lw (LW-1) pt 4 linecolor rgb "black"
set style line 11 lt 5 lw (LW-1) pt 4 linecolor rgb "black"



#parameters from halle
alpha2 = 0.000717635
n2 = 5.5

alpha3 = 0.001223242
n3 = 5.5

rhog=98.1;
set samples 1000
#drainage
pc_vg2(x)=1/alpha2*(x**(-n2/(n2-1))-1)**(1/n2)
#bewasser
pc_vg3(x)=1/alpha3*(x**(-n3/(n3-1))-1)**(1/n3)