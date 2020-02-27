c----------------------------------------------------------------------------------------------
       SUBROUTINE pshk(SS,n,rdn,x1,x2,WPSn,Jacob,IFLG)
       DOUBLE PRECISION P(5,50),M(2,50)
       COMMON /PARTCL/P,M
       REAL*8 PTcut,AYEcut,QCut,Pi
       COMMON /cut/PTcut,AYEcut,QCut,Pi
       DOUBLE PRECISION KPT2(20),KETA(20),KPhi(20),KPT(20)
       COMMON /PARTCL2/ KPT2, KETA, KPhi,KPT
       INTEGER I,J,IFLG
       DOUBLE PRECISION rdn(30),KEPT2(20)
       DOUBLE PRECISION spt1,spt2,E00,PZ0,cutpt,cutpt2,cuteta,cuteta2,cuteta3
       DOUBLE PRECISION SS, sh, x1,x2,pipt2,WPSn,pau2,Jacob



        WPSn   = 0.0D0
        IFLG = 0 
        pipt2=1.0d0
        pau2=(dsqrt(SS)/2.0d0)**2
        cutpt=PTcut
        cutpt2=PTcut
        cuteta=AYEcut
        cuteta2=AYEcut
        cuteta3=AYEcut
        spt1=0.0d0 
        spt2=0.0d0 



       KEPT2(n-1)=dlog(cutpt2**2)+(dlog(pau2)-dlog(cutpt2**2))*rdn(3*n-4)
       KETA(n-1)=-cuteta3+2.0d0*cuteta3*rdn(3*n-5)
       KPT2(n-1)=dexp(KEPT2(n-1))
       KPT(n-1)=dsqrt(KPT2(n-1))
       KPhi(n-1)=2.0d0*PI*rdn(3*n-3)
       spt1=spt1-KPT(n-1)*dcos(KPhi(n-1))
       spt2=spt2-KPT(n-1)*dsin(KPhi(n-1))

       if(n.gt.2) then   
       do I=1,n-2
       KETA(i)=-cuteta+2.0d0*cuteta*rdn(3*i-2)
       KEPT2(i)=dlog(cutpt**2)+(dlog(pau2)-dlog(cutpt**2))*rdn(3*i-1)
       KPT2(i)=dexp(KEPT2(i))
       KPT(i)=dsqrt(KPT2(i))
       KPhi(i)=2.0d0*PI*rdn(3*i)
       spt1=spt1-KPT(i)*dcos(KPhi(i))
       spt2=spt2-KPT(i)*dsin(KPhi(i))
       enddo
       endif

       KETA(n)=-cuteta2+2.0d0*cuteta2*rdn(3*n-2)
       KPT2(n)=spt1**2+spt2**2
       KPT(n)=dsqrt(spt1**2+spt2**2)


           IF(spt1.eq.0.d0) THEN 
                 IF(spt2.GT.0.d0)   KPhi(n) = PI/2.d0
                 IF(spt2.EQ.0.d0)  KPhi(n) = 0.d0
                 IF(spt2.LT.0.d0)    KPhi(n) = 3.*PI/2.d0
           ELSEIF(spt1.LT.0.d0.AND.spt2.GT.0.d0) THEN
                 KPhi(n) = DATAN(spt2/spt1)+PI
           ELSEIF(spt1.LT.0.d0.AND.spt2.LT.0.d0)THEN
                 KPhi(n) = DATAN(spt2/spt1)+PI
           ELSEIF(spt1.GT.0.d0.AND.spt2.LT.0.d0)THEN
                 KPhi(n) = DATAN(spt2/spt1)+2.0d0*PI
           ELSE
                 KPhi(n)= DATAN(spt2/spt1)
           ENDIF
       do I=1,n-1
       pipt2=pipt2*KPT2(i)
       enddo

       E00=0.0d0
       PZ0=0.0d0
       do J=1,n
       P(4,J+2)=dsqrt(KPT2(J)+M(2,J+2))*(dexp(KETA(J))+dexp(-KETA(J)))/2.0d0
       P(3,J+2)=dsqrt(KPT2(J)+M(2,J+2))*(dexp(KETA(J))-dexp(-KETA(J)))/2.0d0
       P(1,J+2)=KPT(J)*dcos(KPhi(j))
       P(2,J+2)=KPT(J)*dsin(KPhi(j)) 
       P(5,J+2)=dsqrt(P(1,J+2)**2+P(2,J+2)**2+P(3,J+2)**2) 
       E00=E00+P(4,J+2)
       PZ0=PZ0+P(3,J+2) 
       enddo

       x1=(E00+PZ0)/dsqrt(SS) 
       x2=(E00-PZ0)/dsqrt(SS) 
       sh=ss*x1*x2
       If(x1.gt.1.0d0.or.x2.gt.1.0d0) then
       IFLG=1
       endif
       P(4,1)=x1*dsqrt(ss)/2.0d0
       P(1,1)=0.d0
       P(2,1)=0.d0
       P(3,1)=x1*dsqrt(ss)/2.0d0
       P(5,1)=x1*dsqrt(ss)/2.0d0
       P(4,2)=x2*dsqrt(ss)/2.0d0
       P(1,2)=0.d0
       P(2,2)=0.d0
       P(3,2)=-x2*dsqrt(ss)/2.0d0
       P(5,2)=x2*dsqrt(ss)/2.0d0

       WPSn=(2.0d0*cuteta)**(n-2)*(2.0d0*cuteta2)*(2.0d0*cuteta3)
     # *(dlog(pau2)-dlog(cutpt**2))**(n-2)*(dlog(pau2)-dlog(cutpt2**2))*(2.*Pi)**(n-1)/4**n*2.0d0
       WPSn=WPSn*0.38937966D+12*(2.0d0*Pi)**(4.0-3.*n)/sh/ss
       Jacob=pipt2

       RETURN                                                                   
       END 

c====================================================
      SUBROUTINE DK2(IY,I1,I2,ZZ,WPS2,Jacob2,IFLG2) 
cc  I2 carries fz I1 carries 1-fz
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
      INTEGER  IFLG2
      DIMENSION zz(2)
      REAL*8 Jacob2
      COMMON /PARTCL/P(5,50),M(2,50)
      REAL*8 PTcut,AYEcut,QCut,Pi
      COMMON /cut/PTcut,AYEcut,QCut,Pi
      REAL*8 PTcut2,AYEcut2,AYEcut3
      COMMON /cut2/PTcut2,AYEcut2,AYEcut3
      REAL*8 lx,lol,uuz,fz,lamd

      WPS2=0.D0 
      IFLG2=0.d0 
      MY=M(1,IY)
      MYSQ=M(2,IY)
      M1SQ=M(2,I1)
      M2SQ=M(2,I2)
      DIFF = MY-M(1,I1)-M(1,I2)
      XLA=((MYSQ-M2SQ-M1SQ)**2-4.*M1SQ*M2SQ)
      lamd=dsqrt(XLA)/MYSQ
      IF(DIFF.LE.0.) GO TO 105
      IF(XLA.GT.0)GO TO 110
105   IFLG2=1 
      RETURN
110   XLA=SQRT(XLA)

      EXCM=(MYSQ+M2SQ-M1SQ)/(2.*MY) 
      PXCM=XLA/(2.*MY)
c.......!!!!!!!!!!!!!!!!!!
      lx=PTcut/P(4,IY)
      IF(lx.GT.1.) then 
      IFLG2=1
      RETURN
      ENDIF  
      lol=dlog(lx/(1.0d0-lx))
      uuz=lol+2*abs(lol)*ZZ(1)
      fz=dexp(uuz)/(1+dexp(uuz))
      CV=(2*fz-1.0d0-M2SQ/MYSQ+M1SQ/MYSQ)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/lamd   
      IF(abs(CV).GT.1.) then 
      IFLG2=1
      RETURN
      ENDIF  

      WPS2=1.0d0*Pi*abs(lol)/dsqrt(1.0d0-MYSQ/P(4,IY)**2)/(2.d0*pi)**(3*2.d0-4.0d0)
      Jacob2=1.0d0/1.0d0/(1.0d0/fz+1.0d0/(1.0d0-fz))
      P(1,I1)=0.
      P(2,I1)=0.
      P(3,I1)=-PXCM 
      P(4,I1)=(MYSQ+M1SQ-M2SQ)/(2.*MY)
      P(1,I2)=0.
      P(2,I2)=0.
      P(3,I2)=PXCM
      P(4,I2)=(MYSQ+M2SQ-M1SQ)/(2.*MY)
      SV=SQRT(ABS(1.-CV*CV))
      CALL ROT12(I1,1,3,SV,CV)
      CALL ROT12(I2,1,3,SV,CV)
      PH = 2.0*PI*ZZ(2) 
      SPH=SIN(PH)
       CPH=COS(PH)
      CALL ROT12(I1,2,1,SPH,CPH)
      CALL ROT12(I2,2,1,SPH,CPH)
C
C     BOOST/ROTATE TO LAB FRAME
C
      ETAY=P(5,IY)/MY
      IF(ETAY.LE.1.E-4)GO TO 200
      GAMMAY=P(4,IY)/MY
      CALL BOOSTZ(I1,GAMMAY,ETAY)
      CALL BOOSTZ(I2,GAMMAY,ETAY)
      CV=P(3,IY)/P(5,IY)
      SV=1.-CV*CV
      IF(SV.GT.-.001)GO TO 130
      WPS2=0. 
      IFLG2=2 
      RETURN
130   SV=SQRT(ABS(SV))
      CALL ROT12(I1,1,3,SV,CV)
      CALL ROT12(I2,1,3,SV,CV)
      PTY=SQRT(P(1,IY)**2+P(2,IY)**2)
      IF(PTY.LE.1.E-4)GO TO 200
      CPH=P(1,IY)/PTY
        SPH=P(2,IY)/PTY
      CALL ROT12(I1,2,1,SPH,CPH)
      CALL ROT12(I2,2,1,SPH,CPH)
200   CALL PSET(I1)
      CALL PSET(I2)
      RETURN
      END

c==============================================================
C-------------------------------THREE-BODY DECAY: AUXILIARY SUBROUTINES
C----------------------------------------------------------------------
C ................. BOOSTN -  general boost ...........................
C
      SUBROUTINE BOOSTN(P,R,Q)
C
C
C     The four vector P is assumed to be given in the rest frame of R,
C     which must be a timelike vector.
C     output Q is the vector P boosted to the frame in which R is given.
C                                              Compare Jackson, p.517
C                                              D. Zeppenfeld (28.6.1985)
C
      REAL*4 P(5),R(4),Q(5)
      REAL*4 BETA(3), X, Y, GAMMA
      INTEGER I

      X = 0D0
      Y = 0D0
      DO I = 1,3
         BETA(I) = R(I)/R(4)
         X = X + BETA(I)**2
         Y = Y + BETA(I)*P(I)
      ENDDO
      IF (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) THEN
         GAMMA = 1D0/DSQRT(1D0-X)
         DO I = 1,3
            Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(4))
         ENDDO
         Q(4) = GAMMA*(P(4) + Y)
      ELSE
         DO I = 1,4
            Q(I) = P(I)
         ENDDO
         IF(X.GE.(1D0-1D-12))
     *      WRITE(6,1000) R,R(4)**2-R(1)**2-R(2)**2-R(3)**2
      ENDIF

      Q(5) = SQRT(Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3))
 1000 FORMAT (' The reference vector ',4G12.3,' is not timelike.'/
     1        ' R**2 = ',G12.3)
      END
C............................................Lorentz boost along Z axis
       SUBROUTINE BOOSTZ(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
         DOUBLE PRECISION M
         COMMON /PARTCL/ P(5,50),M(2,50)

         TMP    = GAMMA*P(3,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(3,I)
         P(3,I) = TMP
         P(5,I) = SQRT(P(1,I)**2 + P(2,I)**2 + P(3,I)**2)
       RETURN 
       END
C............................................Lorentz boost along X axis
       SUBROUTINE BOOST1(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(2,50)

         TMP    = GAMMA*P(1,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(1,I)
         P(1,I) = TMP
         P(5,I) = SQRT(P(1,I)*P(1,I) + P(2,I)*P(2,I) + P(3,I)*P(3,I))
       RETURN
       END
C............................................Lorentz boost along Y axis
       SUBROUTINE BOOSTY(I,GAMMA,ETA)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(2,50)

         TMP    = GAMMA*P(2,I) + ETA*P(4,I)
         P(4,I) = GAMMA*P(4,I) + ETA*P(2,I)
         P(2,I) = TMP
         P(5,I) = SQRT(P(1,I)*P(1,I) + P(2,I)*P(2,I) + P(3,I)*P(3,I))
       RETURN
       END
C......................................................Spatial rotation
       SUBROUTINE ROT12(I,K1,K2,S,C)

         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M 
         COMMON /PARTCL/ P(5,50),M(2,50)

         TMP     = C*P(K1,I) + S*P(K2,I)
         P(K2,I) = C*P(K2,I) - S*P(K1,I)
         P(K1,I) = TMP
       RETURN
       END
C....................................Store magnitude of momentum vector
       SUBROUTINE PSET(I)
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         DOUBLE PRECISION M(2,50),P(5,50),TMP
         COMMON /PARTCL/ P,M

         TMP    = P(1,I)**2 + P(2,I)**2 + P(3,I)**2
         P(5,I) = SQRT(TMP)
       RETURN
       END
C=====================TRANSFER momentums for HELAS===================#

         SUBROUTINE TRANSFER(I,Q)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         COMMON/PARTCL/P(5,50),M(2,50)
         DIMENSION Q(0:3)
                
          Q(0) = P(4,I) 
          Q(1) = P(1,I) 
          Q(2) = P(2,I)
          Q(3) = P(3,I) 

         RETURN
         END
C=====================TRANSFER momentums for HELAS===================#
         SUBROUTINE TRANSFER2(Q,p1,p2,p3,p4,p5)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         DOUBLE PRECISION Q(0:3,5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
           Q(0,1)=p1(0) 
           Q(1,1)=p1(1) 
           Q(2,1)=p1(2) 
           Q(3,1)=p1(3) 
           Q(0,2)=p2(0) 
           Q(1,2)=p2(1) 
           Q(2,2)=p2(2) 
           Q(3,2)=p2(3) 
           Q(0,3)=p3(0) 
           Q(1,3)=p3(1) 
           Q(2,3)=p3(2) 
           Q(3,3)=p3(3) 
           Q(0,4)=p4(0) 
           Q(1,4)=p4(1) 
           Q(2,4)=p4(2) 
           Q(3,4)=p4(3) 
           Q(0,5)=p5(0) 
           Q(1,5)=p5(1) 
           Q(2,5)=p5(2) 
           Q(3,5)=p5(3) 
         RETURN
         END
C=============================Rapidity================================#
      FUNCTION RPID(P)
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
      double precision P(0:3)

      RPID = 0.5d0*dlog((p(0)+p(3))/(p(0)-p(3)))

      return  
      end
C=========================== phi=========================#
         SUBROUTINE COLAT2(I,PH)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         COMMON/PARTCL/P(5,50),M(2,50)
         PI = DACOS(-1.d0)
    
            IF(P(1,I).NE.0.d0) GO TO 88
                 IF(P(2,I).GT.0.d0)  PH = PI/2.d0
                 IF(P(2,I).EQ.0.d0)  PH = 0.d0
                 IF(P(2,I).LT.0.d0)  PH = 3*PI/2.d0
                 RETURN

88            IF(P(1,I).LT.0.d0)THEN
                 PH = DATAN(P(2,I)/P(1,I))+PI
             ELSE
              IF(P(1,I).GT.0.d0.AND.P(2,I).LT.0.d0)THEN
                 PH = DATAN(P(2,I)/P(1,I))+2.*PI
             ELSE
                 PH= DATAN(P(2,I)/P(1,I))
             ENDIF
             ENDIF
             RETURN
             END 
C=========================== phi=========================#
         SUBROUTINE COLAT(P,PH)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         double precision P(0:3)
         PI = DACOS(-1.d0)
    
            IF(P(1).NE.0.d0) GO TO 88
                 IF(P(2).GT.0.d0)  PH = PI/2.d0
                 IF(P(2).EQ.0.d0)  PH = 0.d0
                 IF(P(2).LT.0.d0)  PH = 3*PI/2.d0
                 RETURN

88            IF(P(1).LT.0.d0)THEN
                 PH = DATAN(P(2)/P(1))+PI
             ELSE
              IF(P(1).GT.0.d0.AND.P(2).LT.0.d0)THEN
                 PH = DATAN(P(2)/P(1))+2.*PI
             ELSE
                 PH= DATAN(P(2)/P(1))
             ENDIF
             ENDIF
             RETURN
             END 
C=======================PT=========================#

         FUNCTION PT(P)
         IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)
         double precision P(0:3)

          PT = dsqrt(p(1)**2 + p(2)**2)
            RETURN
            END 

C=============================Rapidity================================#

      SUBROUTINE RAPIDITY(L,AYE)
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z) 
      COMMON/PARTCL/P(5,50),M(2,50)

      aye = 0.5d0*dlog((p(4,l)+p(3,l))/(p(4,l)-p(3,l)))
                                !aye = abs(aye)
      return  
      end
