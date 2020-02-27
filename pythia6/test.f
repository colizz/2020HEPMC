       program test 

c----------------Initialization---------------
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      Double precision MMee
C...Parameters. 
      Common/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)  
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
c***************************************************************************
c      MSTP(41)=0  ! No Decay

c      Z' Decay Channel
      MDME(289,1)=0
      MDME(290,1)=0
      MDME(291,1)=0
      MDME(292,1)=0
      MDME(293,1)=0
      MDME(294,1)=0
      MDME(295,1)=0
      MDME(296,1)=0
      MDME(297,1)=1
      MDME(298,1)=0
      MDME(299,1)=0
      MDME(300,1)=0
      MDME(301,1)=0
      MDME(302,1)=0
      MDME(303,1)=0
      MDME(304,1)=0
      MDME(305,1)=0
      MDME(306,1)=0
      MDME(307,1)=0
      MDME(308,1)=0
      MDME(309,1)=0
      MDME(310,1)=0
c     Z decay channel
      MDME(174,1)=0
      MDME(175,1)=0
      MDME(176,1)=0
      MDME(177,1)=0
      MDME(178,1)=0
      MDME(179,1)=0
      MDME(180,1)=0
      MDME(181,1)=0
      MDME(182,1)=1
      MDME(183,1)=0
      MDME(184,1)=0
      MDME(185,1)=0
      MDME(186,1)=0
      MDME(187,1)=0
      MDME(188,1)=0
      MDME(189,1)=0
c     Gamma decay channel
      MDME(162,1)=0
      MDME(163,1)=0
      MDME(164,1)=0
      MDME(165,1)=0
      MDME(166,1)=0
      MDME(167,1)=0
      MDME(168,1)=0
      MDME(169,1)=0
      MDME(170,1)=1
      MDME(171,1)=0
      MDME(172,1)=0
      MDME(173,1)=0
c*******************************************************************

      MSTP(61)=0  ! No IRS
      MSTP(71)=0  ! No FRS
      MSTP(81)=0  ! No Beam Frag
      MSTP(111)=0 ! No Hadronization
      ECM=14000.D0
      PMAS(32,1)=200.d0
      MSTP(44)=7  ! 7 for full gamma, Z, Z' included
      CKIN(1)=20.d0
      CKIN(2)=-1.0d0
      MSTP(32)=4  ! Choice of Q, Q^2
      
      NEV=10000 

       MSEL=0
       MSUB(141)=1
 
c      CALL PYGIVE('MSEL=21')

c      CALL PYGIVE('MSEL=6')
c      CALL PYGIVE('PMAS(6,1)=175.0')

      CALL PYINIT('CMS','p','p',ECM)

c---------------EVENT LOOP--------------------
      CALL PYBOOK(1,'ME',50,0.D0,250.D0)

      DO 2000 IEV=1,NEV
      CALL PYEVNT
      CALL PYLIST(2)  ! List Event information

c      PT=PYP(11,10)
c      CALL PYFILL(1,PT,1.D0) ! Draw PT figure

       MMee=-(PYP(11,1)+PYP(12,1))**2-(PYP(11,2)+PYP(12,2))**2
     $ -(PYP(11,3)+PYP(12,3))**2+(PYP(11,4)+PYP(12,4))**2

       MMee=dsqrt(MMee)

       CALL PYFILL(1,MMee,1.D0) ! Draw MMee figure

2000  Continue
      
c      CALL PYLIST(12) ! Decay list

c--------------FINALIZATION--------------------
      CALL PYSTAT(1)     

      CALL PYHIST

      print *, PARI(21),PARI(22)

      end
