      SUBROUTINE m2(P,ANS)       
      IMPLICIT NONE
      include "coupl.inc"

      INTEGER                 NCOMB        
      PARAMETER (             NCOMB= 108)
      INTEGER                 NEXTERNAL        
      PARAMETER (             NEXTERNAL=  5)

      REAL*8 P(0:3,NEXTERNAL)
      REAL*8 ANS,T,IDEN,MATRIX

      INTEGER NHEL(NEXTERNAL,NCOMB)

      INTEGER I, IHEL, IHEL2, IC(NEXTERNAL)
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 0/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 0, 0/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1,-1, 1, 0/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 0,-1,-1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 0,-1, 0/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 0,-1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1,-1, 0, 0,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1,-1, 0, 0, 0/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1,-1, 0, 0, 1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1,-1, 0, 1,-1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1,-1, 0, 1, 0/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1,-1, 0, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1,-1, 1,-1, 0/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1,-1, 1, 0, 0/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) /-1,-1, 1, 1, 0/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) /-1, 1,-1,-1, 0/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) /-1, 1,-1, 0, 0/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) /-1, 1,-1, 1, 0/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) /-1, 1, 0,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) /-1, 1, 0,-1, 0/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) /-1, 1, 0,-1, 1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) /-1, 1, 0, 0,-1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) /-1, 1, 0, 0, 0/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) /-1, 1, 0, 0, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) /-1, 1, 0, 1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) /-1, 1, 0, 1, 0/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) /-1, 1, 0, 1, 1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) /-1, 1, 1,-1, 0/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  49),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  50),IHEL=1, 5) /-1, 1, 1, 0, 0/
      DATA (NHEL(IHEL,  51),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  52),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  53),IHEL=1, 5) /-1, 1, 1, 1, 0/
      DATA (NHEL(IHEL,  54),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  55),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  56),IHEL=1, 5) / 1,-1,-1,-1, 0/
      DATA (NHEL(IHEL,  57),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  58),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  59),IHEL=1, 5) / 1,-1,-1, 0, 0/
      DATA (NHEL(IHEL,  60),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  61),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  62),IHEL=1, 5) / 1,-1,-1, 1, 0/
      DATA (NHEL(IHEL,  63),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  64),IHEL=1, 5) / 1,-1, 0,-1,-1/
      DATA (NHEL(IHEL,  65),IHEL=1, 5) / 1,-1, 0,-1, 0/
      DATA (NHEL(IHEL,  66),IHEL=1, 5) / 1,-1, 0,-1, 1/
      DATA (NHEL(IHEL,  67),IHEL=1, 5) / 1,-1, 0, 0,-1/
      DATA (NHEL(IHEL,  68),IHEL=1, 5) / 1,-1, 0, 0, 0/
      DATA (NHEL(IHEL,  69),IHEL=1, 5) / 1,-1, 0, 0, 1/
      DATA (NHEL(IHEL,  70),IHEL=1, 5) / 1,-1, 0, 1,-1/
      DATA (NHEL(IHEL,  71),IHEL=1, 5) / 1,-1, 0, 1, 0/
      DATA (NHEL(IHEL,  72),IHEL=1, 5) / 1,-1, 0, 1, 1/
      DATA (NHEL(IHEL,  73),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  74),IHEL=1, 5) / 1,-1, 1,-1, 0/
      DATA (NHEL(IHEL,  75),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  76),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  77),IHEL=1, 5) / 1,-1, 1, 0, 0/
      DATA (NHEL(IHEL,  78),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  79),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  80),IHEL=1, 5) / 1,-1, 1, 1, 0/
      DATA (NHEL(IHEL,  81),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  82),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  83),IHEL=1, 5) / 1, 1,-1,-1, 0/
      DATA (NHEL(IHEL,  84),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  85),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  86),IHEL=1, 5) / 1, 1,-1, 0, 0/
      DATA (NHEL(IHEL,  87),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  88),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  89),IHEL=1, 5) / 1, 1,-1, 1, 0/
      DATA (NHEL(IHEL,  90),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  91),IHEL=1, 5) / 1, 1, 0,-1,-1/
      DATA (NHEL(IHEL,  92),IHEL=1, 5) / 1, 1, 0,-1, 0/
      DATA (NHEL(IHEL,  93),IHEL=1, 5) / 1, 1, 0,-1, 1/
      DATA (NHEL(IHEL,  94),IHEL=1, 5) / 1, 1, 0, 0,-1/
      DATA (NHEL(IHEL,  95),IHEL=1, 5) / 1, 1, 0, 0, 0/
      DATA (NHEL(IHEL,  96),IHEL=1, 5) / 1, 1, 0, 0, 1/
      DATA (NHEL(IHEL,  97),IHEL=1, 5) / 1, 1, 0, 1,-1/
      DATA (NHEL(IHEL,  98),IHEL=1, 5) / 1, 1, 0, 1, 0/
      DATA (NHEL(IHEL,  99),IHEL=1, 5) / 1, 1, 0, 1, 1/
      DATA (NHEL(IHEL, 100),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL, 101),IHEL=1, 5) / 1, 1, 1,-1, 0/
      DATA (NHEL(IHEL, 102),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL, 103),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL, 104),IHEL=1, 5) / 1, 1, 1, 0, 0/
      DATA (NHEL(IHEL, 105),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL, 106),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL, 107),IHEL=1, 5) / 1, 1, 1, 1, 0/
      DATA (NHEL(IHEL, 108),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA IDEN /216.0d0/


      DO IHEL=1,NEXTERNAL
         IC(IHEL) = +1
      ENDDO

      T=0.0d0       
      DO IHEL2=1,NCOMB
      T=T+MATRIX(P ,NHEL(1,IHEL2),IC(1)) 
      ENDDO
      ANS=T/IDEN
      RETURN
      END
       
      REAL*8 FUNCTION MATRIX(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : uu~ -> ZZZ  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  9,NEIGEN=  1) 
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   17, NCOLOR=  1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      INTEGER                 NEXTERNAL        
      PARAMETER (             NEXTERNAL=  5)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
      include "coupl.inc"
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     3/                                  
C               T[ 2, 1]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL VXXXXX(P(0,3   ),ZMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),ZMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
      CALL FVOXXX(W(1,2   ),W(1,3   ),GZU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL FVOXXX(W(1,6   ),W(1,4   ),GZU ,ZERO    ,ZERO    ,W(1,7   ))    
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GZU ,AMP(1   ))            
      CALL FVIXXX(W(1,1   ),W(1,4   ),GZU ,ZERO    ,ZERO    ,W(1,8   ))    
      CALL FVOXXX(W(1,2   ),W(1,5   ),GZU ,ZERO    ,ZERO    ,W(1,9   ))    
      CALL IOVXXX(W(1,8   ),W(1,9   ),W(1,3   ),GZU ,AMP(2   ))            
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GZU ,AMP(3   ))            
      CALL FVOXXX(W(1,2   ),W(1,4   ),GZU ,ZERO    ,ZERO    ,W(1,10  ))    
      CALL FVIXXX(W(1,1   ),W(1,5   ),GZU ,ZERO    ,ZERO    ,W(1,11  ))    
      CALL IOVXXX(W(1,11  ),W(1,10  ),W(1,3   ),GZU ,AMP(4   ))            
      CALL FVIXXX(W(1,1   ),W(1,3   ),GZU ,ZERO    ,ZERO    ,W(1,12  ))    
      CALL FVIXXX(W(1,12  ),W(1,4   ),GZU ,ZERO    ,ZERO    ,W(1,13  ))    
      CALL IOVXXX(W(1,13  ),W(1,2   ),W(1,5   ),GZU ,AMP(5   ))            
      CALL IOVXXX(W(1,12  ),W(1,10  ),W(1,5   ),GZU ,AMP(6   ))            
      CALL JIOXXX(W(1,1   ),W(1,2   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,14  ))    
      CALL HVVXXX(W(1,4   ),W(1,14  ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     15  ))                                                          
      CALL VVSXXX(W(1,5   ),W(1,3   ),W(1,15  ),GZZH ,AMP(7   ))           
      CALL HVVXXX(W(1,14  ),W(1,3   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     16  ))                                                          
      CALL VVSXXX(W(1,5   ),W(1,4   ),W(1,16  ),GZZH ,AMP(8   ))           
      CALL HVVXXX(W(1,4   ),W(1,3   ),GZZH ,HMASS   ,HWIDTH  ,W(1,         
     &     17  ))                                                          
      CALL VVSXXX(W(1,5   ),W(1,14  ),W(1,17  ),GZZH ,AMP(9   ))           
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)
      MATRIX = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATRIX =MATRIX+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      END
       
       
