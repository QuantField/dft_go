set terminal png                  # Set the size of png
set output "ShrodLog.png" 
set title  'Hydrogen Probability Densities '
 set key left top
#set xlabel 'r'
#set ylabel 'r*psi'
set grid 
set logscale x
# scaling by 0.5 (0.5*$2)
plot  "psi_1.dat" using 1:2 title 'n=1' with lines , \
      "psi_2.dat" using 1:2 title 'n=2' with lines , \
      "psi_3.dat" using 1:2 title 'n=3' with lines 
      
      
