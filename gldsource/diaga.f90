
! diaga.f diagnostics for EMBM atmosphere and sea ice

      subroutine diaga

      include 'var.f90'

      real amin,amax,sum1,sum2,sum3,vsc

      integer i,j,k,iamin,iamax,jamin,jamax

      print*

      call aminmax(imax,jmax,tq(1,1,1),amin,amax,iamin,iamax ,jamin,jamax,2,1)
      print*,'min atm T ',amin,' at ',iamin,jamin
      print*,'max atm T ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,tq(1,1,1),amin,amax,iamin,iamax ,jamin,jamax,2,2)
      print*,'min atm q ',1e3*amin,' at ',iamin,jamin
      print*,'max atm q ',1e3*amax,' at ',iamax,jamax

      call aminmax(imax,jmax,varice,amin,amax,iamin,iamax ,jamin,jamax,2,1)
      print*,'min h ice ',amin,' at ',iamin,jamin
      print*,'max h ice ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,varice,amin,amax,iamin,iamax ,jamin,jamax,2,2)
      print*,'min A ice ',amin,' at ',iamin,jamin
      print*,'max A ice ',amax,' at ',iamax,jamax
             !pptn(maxi,maxj)   ; precipitation
      call aminmax(imax,jmax,pptn(1,1),amin,amax,iamin,iamax ,jamin,jamax,1,1)
      print*,'min pptn  ',amin,' at ',iamin,jamin
      print*,'max pptn  ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,evap(1,1),amin,amax,iamin,iamax ,jamin,jamax,1,1)
      print*,'min evap  ',amin,' at ',iamin,jamin
      print*,'max evap  ',amax,' at ',iamax,jamax

      call aminmax(imax,jmax,pme(1,1),amin,amax,iamin,iamax ,jamin,jamax,1,1)
      print*,'min P-E   ',amin,' at ',iamin,jamin
      print*,'max P-E   ',amax,' at ',iamax,jamax

! compute total water content of planet (should be numerically const.)

      sum1=0.
      sum2=0.
      sum3=0.

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(2,i,j) !specific humidity (2) in atmosphere
            !varice(1,i,j)- sea ice variables: average height (1)
            sum2 = sum2 + varice(1,i,j)
            do k=1,kmax
               sum3 = sum3 + ts(2,i,j,k)*dz(k) !salinity(2)arr
            enddo
         enddo
      enddo

      vsc = ds*dphi*rsc*rsc
! print*,'total water (m^3)',
!hatmbl(2) = 1800. - atmos boundary layer height (m) for heat (1) and humidity (2)
!rhoio - ratio of sea-ice to ocean densities; dsc = 5.e3 - ocean depth scale (m)
!saln0 = 34.9 - reference salinity (can be used as offset)
      print*,'total water', (sum1*rhoao*hatmbl(2) + sum2*rhoio - sum3*dsc/saln0)
! 2     *vsc
      print*,sum1*rhoao*hatmbl(2),sum2*rhoio,-sum3*dsc/saln0

! write(56,*)
! 1    (sum1*rhoao*hatmbl(2) + sum2*rhoio - sum3*dsc/saln0)
! 2, sum2*rhoio

! compute total heat content of planet (should be numerically const.)

      sum1=0.
      sum2=0.

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq(1,i,j)
            do k=1,kmax
               sum2 = sum2 + ts(1,i,j,k)*dz(k)
            enddo
         enddo
      enddo
!ghs - diagnostic global heat source from radiation and
!      water phase change (sb surflux.f)
      vsc = ds*dphi*rsc*rsc*1.e-21
      print*,'total heat (W*10^21)', vsc* (sum1*rhoair*cpa*hatmbl(1) + sum2*dsc*rh0sc*cpsc - ghs*dt(kmax)*tsc) ,sum1*rhoair*cpa*hatmbl(1)*vsc,sum2*dsc*rh0sc*cpsc*vsc ,-ghs*vsc*dt(kmax)*tsc

      end

      subroutine aminmax(imax,jmax,a,amin,amax,iamin,iamax ,jamin,jamax,lmax,l)

      implicit none

      integer i,j,imax,jmax,iamin,iamax,jamin,jamax,lmax,l
      real amin,amax,a(lmax,imax,jmax)

      amin = a(l,1,1)
      amax = a(l,1,1)
      iamin = 1
      iamax = 1
      jamin = 1
      jamax = 1

      do j=1,jmax
         do i=1,imax
            if(a(l,i,j).lt.amin)then
               amin = a(l,i,j)
               iamin = i
               jamin = j
            endif
            if(a(l,i,j).gt.amax)then
               amax = a(l,i,j)
               iamax = i
               jamax = j
            endif
         enddo
      enddo

      end

