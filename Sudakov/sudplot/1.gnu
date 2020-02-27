set pointsize 1.0
set tics in
set mxtics 
set mytics 
set notitle 

set log y
set mytics 10
set terminal  epslatex color
set output "data.tex"
set log x
set title '$Q_0=10$ GeV'

set ylabel '$\Delta(Q^2)$'
set xlabel 'Q (GeV)'
set yrange[0.1:1.1]
plot "herwig"   index 0 u ($1):($2) title 'Herwig' w l  lt 1  ,\
     "nll1"   index 0 u ($1):($2) title 'NLL'  w l lt 2 ,\
     "nll2"   index 0 u ($1):($2) title 'NLL-Appro.'  w l lt 3 

     
set output
!latex tot.tex
!dvips -E tot.dvi -o p2.eps
!epstopdf p2.eps
