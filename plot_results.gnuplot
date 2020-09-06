set title  'Hydrogen (n=1, l=0) radial WF'
set xlabel 'r'
#set ylabel 'r*psi'
set grid 
set logscale x
# scaling by 0.5 (0.5*$2)
plot  "psi.dat" using 1:2 title 'r*psi' with lines, \
      "pot.dat" using 1:(0.01*$2) title 'V(scaled *0.01)' with lines 
      
      
