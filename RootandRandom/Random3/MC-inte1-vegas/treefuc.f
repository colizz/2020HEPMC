      double precision function fes(xx,wgt)
      implicit double precision (a-h,m,o-z)
      dimension xx(10)
      double precision xxx1

      xxx1=xx(1)
      fes=4.d0/(1.d0+xxx1*xxx1)
c      print *, "yy", xx(1), fes
      return
                  
      end
