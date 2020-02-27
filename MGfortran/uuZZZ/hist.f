c============================================================
c
c  A set of subroutines which consecutively initiate a
c  histogram session, add a point to a histogram suite, and
c  end a histogram session.
c
c  Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c  Initial version:  1992 Oct 1  by Adam Duff
c  Rewritten:        1998 Jan 16 by David Rainwater
c
c============================================================
c
c.............................................................
      subroutine hbookdp1(id,name,jbin,xlo,xup,r_zero)
      implicit none
      character*(*) name
      integer id,jbin
      real xlo, xup, r_zero
      include 'histcb.inc'
c
      cfile(id) = name
      
      call hbook1( id,name,jbin,xlo,xup,r_zero)
      end
c ............................................................

      SUBROUTINE XCURVE
      IMPLICIT none 
      integer  nput
c
c declare histogram common variables
c
      include 'histcb.inc'
      integer id, ibin, n, i1, imode
      character kcase*1
      character*72 cnfile(30),name
      real delx,xres(1000,2)
cc
      open(21,file='hist_nos.dat',status='old')

      read(21,*) nput
      if ( nput.gt.30 ) nput = 30
      imode = 0
      do n = 1,nput
         read(21,*) idlist(n)
         id = idlist(n)
      end do

      kcase = '1'
      if ( imode.eq.0 ) then
      NAME = 'dist'
      call selectfile(NAME,'.top')
      OPEN(UNIT=20,FILE=NAME,STATUS='NEW')
      end if

      DO  n = 1,nput
         id = idlist(n)
         print*,' id = ',id,' nbin = ',nbin(id)
         delx = (hup(id)-hlo(id))/(float(nbin(id))+1e-10)
         print*,' delx = ',delx
         XRES(1,1) = hlo(id) + DELX/2.
         DO  IBIN = 2,nbin(id)
            XRES(IBIN,1) = XRES(IBIN-1,1) + DELX
         END DO
         kcase = '1'
         CALL HUNPAK(ID,XRES(1,2),kcase,0)
         do i1 = 1,nbin(id)
            if (xres(i1,2).eq.0.) xres(i1,2) = 1e-10
         end do
         if ( imode.ne.0 ) then
            NAME = CNFILE(n)(1:(index(cnfile(n),' ')-1))//'.top'
            OPEN(UNIT=20,FILE=NAME,STATUS='NEW')
         end if

         write (20,'(A)') ' set intensity 4 '
         write (20,'(A)') ' set font duplex ' 
         write (20,'(A)') ' set tick size .05 '
c         write (20,'(A)') ' set title size 2.0'
         write (20,'(A,i3,a)') ' title 3.0 8.5 ''Fig.',
     1                         id,': '//cfile(id)
         write (20,'(A)') ' set window  x 2.5 to 9.0 y 3.0 to 8.0 '
         write (20,'(A)') '( set scale x lin y log ( lin'
         write (20,'(A)') '( set limits x   to   y   to  '

         DO  IBIN = 1,NBIN(ID)
            WRITE (20,'(2g20.6)')  XRES(ibin,1), xres(ibin,2)
         END DO
         write (20,'(A)') ' hist solid'          

         if ( imode.eq.0 ) then
            if ( n.lt.nput ) then
               write (20,*)
               write (20,'(a)') ' new frame '
            end if
         else
            CLOSE(20)
         end if
      END DO


      if ( imode.eq.0 ) then
         close(20)
      end if
ccc
      return
      end

ccccccccccccccccccccc begin subroutine SelectFile ccccccccccccccc
c
c  This subroutine selects a file name which does not exist yet
c  Terrance Figy <terrance@pheno.physics.wisc.edu
c  August 17, 2001  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine SelectFile(name,extens)
      implicit none
      character*50 name,lname,lname1
      character*1 digit(0:9)
      character*4 extens
      integer ii,i,j
      logical FileAlreadyExists
      data DIGIT /'0','1','2','3','4','5','6','7','8','9'/
 
      lname1 = name
      do 20 ii = 1,99
     
         if((ii.ge.1).and.(ii.le.9)) then
            
            lname = 
     &           lname1(1:index(lname1,'  ')-1)//'.'//DIGIT(ii)//extens
            inquire(file=lname,exist=FileAlreadyExists)
            if(FileAlreadyExists) then
               goto 20
            else
               name = lname(1:index(lname,' ')-1)
               return
            endif
         elseif(ii.gt.9) then
            i = (ii/10)
            j = mod(ii,10)
            lname = lname1(1:index(lname1,'  ')-1)//'.'//
     &              DIGIT(i)//DIGIT(j)//extens
            inquire(file=lname,exist=FileAlreadyExists)
            if(FileAlreadyExists) then
               goto 20
            else
               name = lname(1:index(lname,' ')-1)
               return
            endif
         endif

 20      continue
         write(*,*) '-------- fatal crash in SelectFile ----------'
         stop
         end
c*************************************************************************
