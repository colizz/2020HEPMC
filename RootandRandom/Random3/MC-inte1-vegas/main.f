c  int_0^1 exp(x/2)dx
      include 'treefuc.f'
      program main
      implicit double precision (a-h,m,o-z)
      parameter (pi=3.14159265358979323846264338328)
      common/bveg1/ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      external fes
      open(10, file='output.txt', status='replace')
c*****************define parameters******************************
      read (*,*) ncall
      write (10,*)'ncall=',ncall

      nprn=-1
      xl(1)=0.0d0
      xu(1)=1.0d0
      itmx=3
      call vegas(1,fes,vfes,sd,chi2a)
      write(10,*)vfes, sd, chi2a
      close(10)
		
      end

	  

	
