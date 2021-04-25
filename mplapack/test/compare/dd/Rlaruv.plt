set title "Rlaruv test (0,1) uniform random"
set term png
set output 'Rlaruv.png'

set xlabel "i-th random number"
set ylabel "(0,1)

n=100 
max=0. #max value
min=1. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 # fill style

plot "Rlaruv.txt" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb"green" notitle


