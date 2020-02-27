cli  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
cli   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
cli   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
cli   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
cli   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
cli   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
cli   which is defined as the bottom quark mass, whenever it can be applied.


      subroutine aspdflib(alphas2,SCALE,iord)
      implicit double precision (a-h,o-z)
      double precision NF,iord
cli     DATA XMC/1.43D0/,XMB/4.30D0/,XMT/100.D0/
      DATA XMC/1.3D0/,XMB/4.50D0/,XMT/100.D0/
      DATA ZEROD/0.D0/,PONED/0.001D0/,ONED/1.D0/,TWOD/2.D0/
 
      LO = 1
      if(iord.ne.0) LO = 2
  
cli     TMAS = 180.0d0
      TMAS = 173.0d0 

      ALPHAS2 = ZEROD
      PI=4.0D0*ATAN(ONED)
      B6  = (33.D0-2.D0*6.D0)/PI/12.D0
      BP6 = (153.D0 - 19.D0*6.D0) / PI / TWOD / (33.D0 - 2.D0*6.D0)
      B5  = (33.D0-2.D0*5.D0)/PI/12.D0
      BP5 = (153.D0 - 19.D0*5.D0) / PI / TWOD / (33.D0 - 2.D0*5.D0)
      B4  = (33.D0-2.D0*4.D0)/PI/12.D0
      BP4 = (153.D0 - 19.D0*4.D0) / PI / TWOD / (33.D0 - 2.D0*4.D0)
      B3  = (33.D0-2.D0*3.D0)/PI/12.D0
      BP3 = (153.D0 - 19.D0*3.D0) / PI / TWOD / (33.D0 - 2.D0*3.D0)
      XLC = TWOD * LOG( XMC/QCDL5)
      XLB = TWOD * LOG( XMB/QCDL5)
      XLT = TWOD * LOG( XMT/QCDL5 * TMAS/XMT)
      XLLC = LOG( XLC)
      XLLB = LOG( XLB)
      XLLT = LOG( XLT)
      C65  =  ONED/( ONED/(B5 * XLT) - XLLT*BP5/(B5 * XLT)**2 ) - 
     &     ONED/( ONED/(B6 * XLT) - XLLT*BP6/(B6 * XLT)**2 )
      C45  =  ONED/( ONED/(B5 * XLB) - XLLB*BP5/(B5 * XLB)**2 ) - 
     &     ONED/( ONED/(B4 * XLB) - XLLB*BP4/(B4 * XLB)**2 )
      C35  =  ONED/( ONED/(B4 * XLC) - XLLC*BP4/(B4 * XLC)**2 ) - 
     &     ONED/( ONED/(B3 * XLC) - XLLC*BP3/(B3 * XLC)**2 ) + C45
  
      Q = SCALE
      IF  ( Q .GT. XMB ) THEN
        NF = 5.D0
        QCDL5=0.165d0
        if(iord.ne.0)   QCDL5=0.226d0
      ELSE
        NF = 4.D0
        QCDL5=0.215d0 
        if(iord.ne.0)   QCDL5=0.326d0
      ENDIF

 

      XLQ = TWOD *  LOG( Q/QCDL5 )
      XLLQ =  LOG( XLQ )

      IF      ( NF .EQ. 6.D0 ) THEN
        ALF = ONED/(ONED/(ONED/(B6*XLQ)- BP6/(B6*XLQ)**2*XLLQ) + C65)
        IF (LO.EQ.1) ALF = ONED/B6/XLQ
      ELSEIF  ( NF .EQ. 5.D0 ) THEN
        ALF = ONED/(B5 * XLQ) -  BP5/(B5 * XLQ)**2 * XLLQ
        IF (LO.EQ.1) ALF = ONED/B5/XLQ
      ELSEIF  ( NF .EQ. 4.D0 ) THEN
        ALF = ONED/(ONED/(ONED/(B4*XLQ)- BP4/(B4*XLQ)**2*XLLQ) + C45)
        IF (LO.EQ.1) ALF = ONED/B4/XLQ
      ELSEIF  ( NF .EQ. 3.D0 ) THEN
        ALF = ONED/(ONED/(ONED/(B3*XLQ)- BP3/(B3*XLQ)**2*XLLQ) + C35)
        IF (LO.EQ.1) ALF = ONED/B3/XLQ
      ELSE
        WRITE(*,*) 'Error in Alphas2'
        STOP
      ENDIF
      ALPHAS2 = ALF
      RETURN
      END 

