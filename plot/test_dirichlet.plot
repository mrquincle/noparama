# Plot file for data coming out of a Dirichlet distribution
# This should just be a bunch of probability distribution. The Dirichlet means that it's a mixture, not much more.
#   @param filename

# Print to stdout
set print "-"

if (!exists("filename")) print "A filename required as argument"

# Number of intervals
n=50 
max=10.
min=0. 

# Width af an interval
width=(max-min)/n 

# Histogram
hist(x,width) = width*(floor((x-min)/width)+0.5) + min

# Graphics
set boxwidth width*0.9

plot filename u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle

pause -1
