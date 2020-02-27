c===================================================================
c
c   Defines common block /histcb/
c
c   Adam Duff & Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c   Initial version:  1992 Oct 1
c   Last modified:    1993 Apr 12
c
c===================================================================
c
      integer  nhmax
      parameter( nhmax=256 )

      integer  nbin(nhmax)
      real*4 hlo(nhmax), hup(nhmax), rwid(nhmax), dbin(nhmax)
      common /histcb/ nbin, hlo, hup, rwid, dbin

      character*72 cfile
      integer idlist
      common /bxcurve/ idlist(30),cfile(nhmax)

