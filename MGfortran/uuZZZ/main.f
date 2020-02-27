      Program main       
   
      IMPLICIT NONE
      include "coupl.inc"
      REAL*8 ANSa,dsigma  

      INTEGER                 NEXTERNAL        
      PARAMETER (             NEXTERNAL=  5)

      REAL*8 ENERGY,SP,S
      REAL*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      REAL*8 Q(0:3,5),Qb(0:3,5)
      REAL*8 X1,X2
       
      REAL*8 P(5,50),M(2,50)
      COMMON /PARTCL/P,M
      
      REAL*8 PT3,PT4,AYE3,AYE4 
      REAL*8 PTcut,AYEcut,QCut,Pi
      COMMON /cut/PTcut,AYEcut,QCut,Pi
      DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
      COMMON /PARTCL2/ KPT2, KETA, KPhi, KPT


      REAL*8 Myrand(30),WGH,hwid
      REAL*8 WPS3,Jacob3
      Integer I,J,IFLAG,IFL,h
      character name*64   
      real*8 f1(-6:6),f2(-6:6),muf,mur
      real*8 alphasPDF


      integer iteration, seed,N_iterations1,N_points,PS_dimension
      integer*8 Ncall
      character*1 digit(9)
      data digit /'1','2','3','4','5','6','7','8','9'/
      character*50 file_name_out, file_name_in,gridname1
      integer max_PS_dim
      parameter (max_PS_dim = 30)
      real*8 rand(max_PS_dim),rn
      real*8 weight,dsig0
      double precision xsection, sdev,chi2
      include 'hist.inc'


C ----------
C BEGIN CODE
C ----------
         gridname1="test"
         file_name_in  = gridname1(1:INDEX(gridname1,'  ')-1)

         pi= DACOS(-1.d0)
         N_iterations1=6.d0
         N_points=20.d0
         PS_dimension=7.d0
  
         ENERGY=14000.0d0
         SP = ENERGY**2
         zmass= 9.11876000e+01 
         zwidth=2.44140351E+00
         hmass=120.0d0
         hwidth=5.75308848E-03
c------------test----------------------
         PTcut=0.001d0
         AYEcut=dlog(2.0d0*ENERGY/zmass)
         print *, AYEcut 
c---------------------------------
         print *, "0000"
c         Call HLIMIT (10000)
         print *, "0001"

c         Call HBOOKdp1(1 ,'d[s]/dPT_1  ',nbin(1 ),hlo(1 ),hup(1 ),0.0)
c         Call HBOOKdp1(2 ,'d[s]/dET_1  ',nbin(2 ),hlo(2 ),hup(2 ),0.0)
         do  I = 1,2
           call hidopt( I, 'inte' )              ! for integrate data
         enddo

        do i = 1, nhmax
         hwid = hup(i) - hlo(i)
         if ( nbin(i).gt.0 .and. hwid.ne.0.0d0 ) then
            rwid(i) = real( nbin(i) ) / hwid
         else
            rwid(i) = 0.0d0
         end if
         end do
         print *, "1111"
         name='NNPDF23_lo_as_0130_qed'
         call InitPDFsetByName(name)
         call InitPDF(0)
          print *, "2222"

         DO iteration = 1, N_iterations1
            Ncall = 2**( N_points - N_iterations1 + iteration )
           
            file_name_out = gridname1(1:INDEX(gridname1,'  ')-1)//
     &                   '.out.'//digit(iteration)
         
          if ( iteration .eq. 1 ) then
               call monaco_init( PS_dimension, Ncall, seed )
               call monaco_read(file_name_in)
            else
               call monaco_init2( PS_dimension, Ncall )
            end if       

          if (iteration.eq.N_iterations1) then
               call monran_set(1) ! save random number info
            endif
 
c************************************************************************
       do J= 1, Ncall

       IF(mod(J,40000).eq.0) then
       print *, J
       endif  

       call monaco_get( rand, weight)
       DO I = 1,30
       Myrand(I) = rand(I)
       ENDDO

 
         bmass=4.2d0
         alpha=1.0d0/1.32506980E+02  
         M(1,1) = 0.d0
         M(1,2) = 0.d0
         M(1,3) = zmass
         M(1,4) = zmass
         M(1,5) = zmass

         DO 50 I = 1,50
         M(2,I) = M(1,I)**2
 50      CONTINUE

c.--------------------------------------------------------------------------------------------------
c....1,2 are the initial partons, 3,4,.. n+1,n+2 are the final states
         CALL  pshk(SP,3,myrand,X1,X2,WPS3,Jacob3,IFLAG)
c.--------------------------------------------------------------------------------------------------

         IF( IFLAG .ne. 0) then
         dsig0 = 0.d0
         goto 1000
         endif

         S = X1*X2*SP

           CALL TRANSFER (1, p1)     
           CALL TRANSFER (2, p2)
           CALL TRANSFER (3, p3)
           CALL TRANSFER (4, p4)
           CALL TRANSFER (5, p5)

           CALL TRANSFER2(Q,p1,p2,p3,p4,p5)

           PT3=KPT(1)
           PT4=KPT(2)
           AYE3=KETA(1)
           AYE4=KETA(2)

      mur=zmass 
      muf=zmass 
      ALFAS=alphasPDF(mur)
      G = DSQRT(4d0*pi*ALFAS) 
      GG(1) = -G
      GG(2) = -G
      Call setpara

      call evolvePDF(X1,muf,f1) 
      call evolvePDF(X2,muf,f2) 

       Call m2(Q,ANSa)

      dsigma= ANSa*f1(2)*f2(-2)/X1/X2
      dsig0=dsigma*WPS3*Jacob3

      if (iteration.eq.N_iterations1) then 
      WGH=weight*dsig0 
c      CALL HF1(1, Sngl(PT3), sngl(WGH)*rwid(1 )) 
c      CALL HF1(2, Sngl(AYE3), sngl(WGH)*rwid(2 )) 
      endif

1000      call monaco_put(rand, weight, dsig0)
         
      enddo 
c***********************************************************************************
      call monaco_write( file_name_out )

      END DO

        call monaco_result(xsection, sdev, chi2)
     
c     print out the result. 
         print *,""
         print *," result (LO): ",xsection," +- ",sdev," fb  "
     &        ,sdev/xsection*100d0,'%'
         print *,""

       Call HISTDO         
       CALL XCURVE 
      end  
