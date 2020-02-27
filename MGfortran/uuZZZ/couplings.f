      subroutine setpara
c***********************************************************************
c This subroutine sets up the HELAS couplings of the MODEL.
c***********************************************************************
      implicit none
      include 'coupl.inc'

      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2,Sqrt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Sqrt2   = 1.414213562d0 )
      double precision  Pi, Fourpi,sin2w,sw,cw,wm,gfermi
      double precision  v
      double precision ez, ey, sc2
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
 
      gfermi= 1.16639000e-05
      ee2 = alpha * Fourpi
      ee  = sqrt( ee2 )
      awidth=0.0d0
      wm = sqrt(zmass**2/2d0+
     $     sqrt(zmass**4/4d0-Pi/Sqrt2*alpha/gfermi*zmass**2))
      sin2w  = 1d0-(wm/zmass)**2
      cw  = sqrt( 1d0 - sin2w )
      sw  = sqrt( sin2w )
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee  
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)

      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw

      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )

      gwwh  = dcmplx( ee2/sin2w*Half*v, Zero )
      gzzh  = dcmplx( ee2/sc2*Half*v, Zero )
      ghhh  = dcmplx( -hmass**2/v*Three, Zero )
      gwwhh = dcmplx( ee2/sin2w*Half, Zero )
      gzzhh = dcmplx( ee2/sc2*Half, Zero)
      ghhhh = ghhh/v

      return
      end

      
