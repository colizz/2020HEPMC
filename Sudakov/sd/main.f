      Program main       
      Implicit double precision (a-h,m,o-z)
      COMMON/BVEG1/NCALL,ITMX,NPRN,NDEV,XL(10),XU(10),ACC
      double precision tcutoff,tup,iord,ischeme
      COMMON/user1/ tcutoff,tup,iord,ischeme
      double precision stcutoff,stup
      double precision sdexp,fsudkov

      external fes
C ----------
C BEGIN CODE
C ----------
        read *,stup
        tup=stup**2
        stcutoff=10.0d0
        tcutoff=stcutoff**2
c....... as LO (iord=0) or NLO (iord ne 0) running
        iord=1.  
c--------CKKW NLL1 (ischeme=0) or Herwig (ischeme=1) or NLL2 (ischeme=2), see hep-ph/0312274
        ischeme=2.
c----------------

	nprn=-1
	ncall=20000
	itmx=4
	xl(1)=0.0d0
	xu(1)=1.0d0
	xl(2)=0.0d0
	xu(2)=1.0d0

	call vegas(2,fes,vfes,sd,chi2a)
	print *,vfes,sd,chi2a
        sdexp=vfes   
        fsudkov=dexp(-sdexp)
	print *, dsqrt(tup), fsudkov
	end

c***************************************************
c***************************************************
	double precision function fes(xx)
        double precision xx(2)
        double precision tcutoff,tup,ischeme
        COMMON/user1/ tcutoff,tup,iord,ischeme
    	double precision it1,it2,zup,zlow,jacob1,jacob2
        double precision itt,zzz,astrong,muscale,iord
        double precision fsplqq,pi
        pi=dacos(-1.0d0)

        if(ischeme.eq.0.) then
        it1=4.0d0*tcutoff
        it2=tup
        itt=it1+(it2-it1)*xx(1)
        jacob1=(it2-it1)
        zlow=dsqrt(tcutoff/itt)
        zup=1.0d0-dsqrt(tcutoff/itt)
        zzz=zlow+(zup-zlow)*xx(2)
        jacob2=(zup-zlow)
        muscale=dsqrt(itt)*zzz*(1.0d0-zzz)
        elseif(ischeme.eq.1.) then
        it1=4.0d0*tcutoff
        it2=tup
        itt=it1+(it2-it1)*xx(1)
        jacob1=(it2-it1)
        zlow=dsqrt(itt/4.0d0/tup)
        zup=1.0d0-dsqrt(itt/4.0d0/tup)
        zzz=zlow+(zup-zlow)*xx(2)
        jacob2=(zup-zlow)
        muscale=dsqrt(itt)
        else
        it1=2.0d0*dsqrt(tcutoff)
        it2=dsqrt(tup)
        itt=it1+(it2-it1)*xx(1)
        jacob1=(it2-it1)
        zlow=0.0d0
        zup=1.0d0
        zzz=zlow+(zup-zlow)*xx(2)
        jacob2=(zup-zlow)
        muscale=itt

        ENDIF


        fsplqq=4.d0/3.d0*(1.d0+zzz**2)/(1.d0-zzz)
        call aspdflib(astrong,muscale,iord)
	  fes=astrong/2.0d0/pi*fsplqq/itt*jacob1*jacob2

        if(ischeme.eq.2) then
        fes=astrong*8.d0/3.d0/pi/muscale*(dlog(it2/muscale)-0.75d0)*jacob1*jacob2
        endif

        end  

