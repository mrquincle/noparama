n=50 #number of intervals
max=10. #max value
min=0. #min value
width=(max-min)/n #interval width
hist(x,width) = width*(floor((x-min)/width)+0.5) + min
set boxwidth width*0.9

plot "test_dirichlet.data" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle

pause -1
