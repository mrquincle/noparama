# Print to stdout
set print "-"

if (!exists("filename")) print "A filename required as argument"

plot filename using 1:2

pause -1
