
! subroutine gseta, sets up  atmosphere and sea-ice

      subroutine gseta

      include 'var.f90'

      real tv, tv2, tv3, tatm, relh0_ocean, relh0_land, circum, sigma
      real scl_co2, u_tau_ice,ch_ice,cpo_ice !pc_co2_rise - excluded ###
      real asurf(maxj), aproj(maxj), obl(maxj)
      real uatm_in(2,maxi,maxj)
      real extra1a,extra1b,extra1c

      integer i, j, k, l, natl1a, npac1a, natl1b, npac1b, natl1c, npac1c
      integer nboxmin, nboxtot, nboxe, nboxw, nboxn, nboxs
	!######### for seasonal run
      integer ii, nn1,nn2
	real SCOSZ(maxj)
	logical end_of_day
	!#########
! constants used in SSW parameterization

      sigma = 5.67e-8
      emo = 0.94 * sigma
      ema = 0.85 * sigma
      solconst = 1368. !W/m**2
      tfreez = 0.
      alpha = 0.3
      alphice = 0.7

! constants used in OLW parameterization

      b00 = 2.43414e2
      b10 = -3.47968e1
      b20 = 1.02790e1
      b01 = 2.60065
      b11 = -1.62064
      b21 = 6.34856e-1
      b02 = 4.40272e-3
      b12 = -2.26092e-2
      b22 = 1.12265e-2
      b03 = -2.05237e-5
      b13 = -9.67000e-5
      b23 = 5.62925e-5

! EMBM stuff follows...
            !наклон
! areas and obliquity factor needed for embm
      ! open(321,file=trim(path_source)//'solar.dat')  !###########

      do j = 1,jmax    !asurf(j)- grid cell area (m**2)
         asurf(j) = rsc*rsc*ds*dphi  !rsc = 6.37e6 - Earth rad (m)
         aproj(j) = rsc*rsc*c(j)*ds*dphi/pi ! c(j)-cos of latitude
         
         obl(j) = 0.25*pi*(1.000-(0.477*((3.*s(j)**2)-1.)/2.))/c(j)
	!solar forcing data array:
         solfor(j) = solconst*obl(j)*aproj(j)/asurf(j)
! print*, j,asurf(j),aproj(j),obl(j), solfor(j)
!#####         write(321,321) j,asurf(j),aproj(j),obl(j), solfor(j)
      enddo
 321  format(1x, i3,4e12.3)
 !############# seasonal solar
!         UPDATE DAY, EARTH-SUN DISTANCE, SOLAR DECLINATION
      end_of_day= (MOD(nint(t*tsc/86400),24).EQ.0)

      IF (end_of_day) then
	 ! sdedy - day of year
	 ! SDEYR - year N
       call SDET_OC(SCOSZ) !mean day insolation W/m**2
        print*, ' time= ',t*tsc/86400
        print*, 'j  solfor(j) SCOSZ(j)'
        write(6,'(i3,2f8.2)')  (j, solfor(j),SCOSZ(j),j=1,jmax)
	! pause 'insolation press enter'
!	 solfor=SCOSZ
	endif
 !#############
! climatological albedo (similar to Weaver et al. 2001)
! or read interpolated data (rma, 20/9/02)

! open(10,file=trim(path_source)//'albedo.dat')
      open(10,file=trim(path_source)//'albedo_function.dat')

      do j=1,jmax
         tv = asin(s(j))  ! latitude (rad)
         tv2 = 0.2 + 0.36*0.5*(1.0 - cos(2.0*tv))
         write(10,*) 180/pi*tv,tv2 ! latitude (grad), albedo(j)
! read(10,*) tv,tv2
         do i=1,imax
            albcl(i,j) = tv2  !albedo(j)
         enddo
      enddo

      close(10)

! atmospheric SSW absorption coefficient, value over land purely diagnostic

      do j=1,jmax
         do i=1,imax
            if(k1(i,j).le.kmax)then
               ca(i,j)=0.3
            else
               ca(i,j)=1.0
            endif
         enddo
      enddo

! read some scalings

      print*,'scl_co2,coef_co2,Time_co2' !###
      read(5,*)scl_co2,coef_co2,Time_co2
!###      scl_co2=1
      print*,scl_co2,coef_co2,Time_co2 !###

! factor corresponding to radiative forcing of 4 W/m**2
! per doubling of atmospheric CO2

      delf2x = 5.77

! present day CO2 concentration (ppm)

      co20 = 350. !1980 year
 !     co20 = 390. !2010 year
      !######CO2 up to 2100
      !B1_CO2(i),A2_CO2(i)- two different CO2 scenario - use one
      open (11, file=trim(path_source)//'CO2B1A2.txt')
       print *, ' Year B1_CO2 A2_CO2'
      do i=1,13
       read (11,*) Year2100(i),B1_CO2(i),A2_CO2(i)
       print '(I6, 2F7.1)', Year2100(i),B1_CO2(i),A2_CO2(i)
      enddo
      istepCO2=0
            
      ! pause 'read CO2 data' !GLD GGGGGGGGG
      
      !######CO2
	if (Time_co2.ne.0) then !####
		delta_co2 = (coef_co2-1)*co20/(Time_co2*365.)
	else
	    delta_co2 = 0.
	end if     !####
! scale co20 by factor scl_co2

      do j = 1,jmax
         do i = 1,imax
            co2(i,j) = scl_co2*co20
         enddo
      enddo
      co2=350 !2100  !<- 1980   A2_CO2(13)

! convert compound annual change in co2 to mks rate of change per dt
!####      ryear = 1./(365.*86400.)
! rate_co2 = (pc_co2_rise/100)*tsc*dtatm*ndta*ryear
!####      rate_co2 = (1. + 0.01*pc_co2_rise)**(tsc*dtatm*ndta*ryear) - 1.

! more constants

      rhoair = 1.25 !air density
      rho0 = 1.e3  !water density (kg/m**3)
      rhoao = rhoair/rho0
! depth scale for atmospheric thermal b.l. used by Weaver et al. (2001)
      hatmbl(1) = 8400. !atmospheric boundary layer height for heat (1)
                       	! and humidity (2)
      cpa = 1004.  !specific heat capacity of air (J/kg/К)
! latent heat of vapourization (J/kg)
      hlv = 2.501e6
! latent heat of sublimation (J/kg)
      hls = 2.835e6
! latent heat of fusion of ice (J/kg)
      hlf = 3.34e5

! scaling for heat forcing of atmosphere
              !usc - velocity scale
      rfluxsca = rsc/(hatmbl(1)*usc*rhoair*cpa)

! atmospheric winds

!c!      open(35,file=trim(path_source)//'uncep.silo')
!c!      read(35,*)((uatm_in(1,i,j),i=1,imax),j=1,jmax)
!c!      close(35)

!c!      open(35,file=trim(path_source)//'vncep.silo')
!c!      read(35,*)((uatm_in(2,i,j),i=1,imax),j=1,jmax)
!c!      close(35)

! conditional zonal average

!c!      do j=1,jmax
!c!         if(j.le.2.or.j.ge.jmax-1)then
!c!         do l=1,2
!c!            tv = 0.
!c!            do i=1,imax
!c!               tv = tv + uatm_in(l,i,j)
!c!            enddo
!c!            tv = tv / imax
!c!            do i=1,imax
!c!               uatm_in(l,i,j) = tv
!c!            enddo
!c!         enddo
!c!         endif
!c!      enddo

! rotate 2 gridboxes to west

!c!      do j=1,jmax
!c!      do l=1,2
!c!       ii=0
!c!       do i=3,imax
!c!         ii=ii+1
!c!         uatm(l,ii,j) = uatm_in(l,i,j)
!c!       enddo
!c!       do i=1,2
!c!        ii=ii+1
!c!        uatm(l,ii,j) = uatm_in(l,i,j)
!c!       enddo
!c!      enddo !l
!c!      enddo !j

!c!      open(35,file=trim(path_source)//'uncep_rot.silo')
!c!	!uatm(2,maxi,maxj) - prescribed atmospheric advective velocities
!c!      write(35,*)((uatm(1,i,j),i=1,imax),j=1,jmax)
!c!      close(35)

!c!      open(35,file=trim(path_source)//'vncep_rot.silo')
!c!      write(35,*)((uatm(2,i,j),i=1,imax),j=1,jmax)
!c!      close(35)

! remove zonal average of v else fail mass conservation (may not be
! disastrous).

!c!      do i=1,imax
!c!         do j=1,jmax                 !usc = 0.05   vel (m/s)
!c!            uatm(1,i,j) = uatm(1,i,j)/usc
!c!            uatm(2,i,j) = uatm(2,i,j)/usc
!c!           uatm(2,i,j) = uatm(2,i,j)*0.0
!c!            uatm(3,i,j) = 0. !not used
!c!         enddo
!c!         uatm(2,i,jmax) = 0.
!c!      enddo

! parameter beta relates vertically-averaged advective transport
! to surface advective transport

! heat advection
      beta(1) = 0.0
! beta(1) = 0.4

! moisture advection
! beta(2) = 0.0
      beta(2) = 0.4

! parameters for extra heat diffusion where pptn high
!diffmod0 - extra heat diffusivity for strong precipitation, set 0 by default
! diffmod0 = 60e6
      diffmod0 = 0.
! ppmin,ppmax  - parameters for above
      ppmin = 2./(365.*86400.)
      ppmax = 4./(365.*86400.)

! nre simpler diffusivity
! or Weaver's data (rma 19/9/02)

! open(47,file=trim(path_source)//'diff.dat')
! open(47,file=trim(path_source)//'diff_DIF.dat')
      open(47,file=trim(path_source)//'diff_ADVDIF.dat')
      do j=1,jmax
         tv = asin(s(j))  !lat
         tv2 = asin(sv(j))
! diffa(1,1,j) = 2.e6
! diffa(1,2,j) = 2.e6
! diffa(1,1,j) = (1.e6 + 2.0*1.e6*0.5*(1. + cos(2.*tv)))
! diffa(1,2,j) = (1.e6 + 2.0*1.e6*0.5*(1. + cos(2.*tv2)))
! diffa(2,1,j) = diffa(1,1,j)
! diffa(2,2,j) = diffa(1,2,j)
! diffa(1,1,j) = (1.e5 + 3.9*1.e6*0.5*(1. + cos(2.*tv)))
! diffa(1,2,j) = (1.e5 + 3.9*1.e6*0.5*(1. + cos(2.*tv2)))
! diffa(2,1,j) = 1.e6
! diffa(2,2,maxj)  - atmospheric diffusivity (l,m,j) where l
!denotes T or Q, m denotes u or v points, for each j point
! diffa(2,2,j) = 1.e6
! Weaver diffusivities as interpolated by Bob
        read(47,'(4x,5e15.5)')tv3,diffa(1,1,j) ,diffa(2,1,j),diffa(1,2,j),diffa(2,2,j)
         diffa(1,1,j) = diffa(1,1,j)/(rsc*usc)
         diffa(1,2,j) = diffa(1,2,j)/(rsc*usc)
         diffa(2,1,j) = diffa(2,1,j)/(rsc*usc)
         diffa(2,2,j) = diffa(2,2,j)/(rsc*usc)
! write(47,'(i4,5e15.5)')j,asin(s(j)),diffa(1,1,j)*rsc*usc
! 1   ,diffa(2,1,j)*rsc*usc,diffa(1,2,j)*rsc*usc,diffa(2,2,j)*rsc*usc
      enddo
      close(47)

! scale height for specific humidity (Peixoto and Oort 1992)
      hatmbl(2) = 1800.

! consts for saturation specific humidity (Bolton 1980)

      const1 = 3.80
      const2 = 21.87
      const3 = 265.5
      const4 = 17.67
      const5 = 243.5

! threshold relative humidity

      rmax = 0.85

! scaling for P-E forcing of atmosphere

! rpmesca = rsc*rho0/(dsc*usc*rhoair)!rsc = 6.37e6 Earth rad (m)
      rpmesca = rsc*rho0/(hatmbl(2)*usc*rhoair)

! reconstruct surface wind field for bulk turbulent transfer and
! zonally average near poles as for uatm for stability

! open(55,file=trim(path_source)//'usurf.dat')

      cd = 0.0013 !drag coefficient for wind stress calculation

      do j=1,jmax
         tv3 = 0.
         do i=1,imax
            if(i.eq.1) then
               tv = (tau(1,i,j)+tau(1,imax,j))/2.
            else
               tv = (tau(1,i,j)+tau(1,i-1,j))/2.
            endif
            if(j.eq.1) then
               tv2 = tau(2,i,j)/2.
            else
               tv2 = (tau(2,i,j)+tau(2,i,j-1))/2.
            endif
            usurf(i,j) = sqrt((sqrt(tv**2 + tv2**2)) *rh0sc*dsc*usc*fsc/(rhoair*cd*scf))
! 1          *rh0sc*dsc*usc*fsc/(rhoair*cd))
            tv3 = tv3 + usurf(i,j)
         enddo
         do i=1,imax
            if(j.le.2.or.j.ge.jmax-1)usurf(i,j) = tv3/imax
! write(55,*)usurf(i,j)
         enddo
      enddo

! close(55)
!c-----------------------------------------------------------------------
! sea ice parameters
!c-----------------------------------------------------------------------

! freezing temperature for average seawater (deg C)
      tsic = -1.8
! constant ice conductivity (W/m/K)
      difsic = 2.166
! in parameterization of heat flux at base of sea ice:
! empirical constant
      ch_ice = 0.0058
! skin friction velocity (m/s)
      u_tau_ice = 0.02
! specific heat of sea water under ice at constant pressure (J/kg/K)
      cpo_ice = 4044.
! representative ice density (kg/m**3)
      rhoice = 913.
! representative under-seaice water density (kg/m**3)
! rho0sea = 1035.
! demarcation thickness between thick and thin ice (m)
      h0sic = 0.01
! useful constant proportional to inverse timscale for surface freezing
! rsictscsf = ch_ice*u_tau_ice*rho0sea*cpo_ice
      rsictscsf = ch_ice*u_tau_ice*rho0*cpo_ice
      print*,'rsictscsf = ',rsictscsf  !dsc = 5e3 depth (m)
      rsictscsf = dsc*dz(kmax)*rho0*cpo_ice/(9.*86400.0)
!      rsictscsf = dsc*dz(kmax)*rho0*cpo_ice/(17.5*86400.0)
      print*,'rsictscsf = ',rsictscsf
! minimum average sea-ice thickness over a grid cell
      hmin = 0.01
      rhmin = 1.0/hmin
! density ratios
      rhooi = rho0/rhoice
      rhoio = rhoice/rho0
! melting factor
      rrholf = 1.0/(rhoice*hlf)

! read initial atmos state
      print*,'tatm relh0_ocean relh0_land'
      read(5,*)tatm,relh0_ocean,relh0_land
!###      tatm=0.
!###	relh0_ocean=0.
!###	relh0_land=0.
      print*,tatm,relh0_ocean,relh0_land

! read freshwater flux perturbation data
      print*,'extra0 range0 nsteps_extra0'
      read(5,*)extra0,range0,nsteps_extra0
!###      extra0=1 !for melt
!###	range0=0.
!###	nsteps_extra0=6250 !for melt
      print*,extra0,range0,nsteps_extra0

! denominator for converting extra0 to mks runoff units:
! over 18 gridboxes within convective regions (50-70N) in N.Atlantic:
      rextra0 = 1.e6/(18.*asurf(1)) !asurf - grid cell area m2
! over 62 gridboxes south of convective regions (20-50N) in N.Atlantic:
!      rextra0 = 1e6/(62.*asurf(1))

! implicit Atlantic-to-Pacific freshwater fluxes in south Atlantic,
! tropical Atlantic and north Atlantic: extra1a, extra1b, extra1c (Sv)
      extra1a = -0.03
      extra1b = 0.17
      extra1c = 0.18

! read scaling factor for extra1a, extra1b, extra1c
      print*,'scl_fwf'
      read(5,*)scl_fwf
!###      scl_fwf=0.75
      print*,scl_fwf

! apply scl_fwf
      extra1a = scl_fwf*extra1a
      extra1b = scl_fwf*extra1b
      extra1c = scl_fwf*extra1c

! use extra1a, extra1b, extra1c, basins data to set up P-E adjustments

! find total no. of Pac/Atl gridboxes

! in south Atlantic (to 20 deg S)
      npac1a = 0 !Pac
      natl1a = 0 !Atl
      do j=16,25
         npac1a = npac1a + ipf(j) - ips(j) + 1
         natl1a = natl1a + iaf(j) - ias(j) + 1
      enddo

! in tropical Atlantic (20 deg S to 24 deg N)
      npac1b = 0
      natl1b = 0
      do j=26,51
         npac1b = npac1b + ipf(j) - ips(j) + 1
         natl1b = natl1b + iaf(j) - ias(j) + 1
      enddo

! in north Atlantic (north of 24 deg N) NB INCLUDES DRY POINTS
      npac1c = 0
      natl1c = 0
      do j=52,jmax
         do i=ips(j),ipf(j)
            if(k1(i,j).le.kmax)npac1c = npac1c + 1
         enddo
         do i=ias(j),iaf(j)
            if(k1(i,j).le.kmax)natl1c = natl1c + 1
         enddo
      enddo

      print*,natl1a, npac1a, natl1b, npac1b, natl1c, npac1c
      print*,'natl1a, npac1a, natl1b, npac1b, natl1c, npac1c '

! increase/decrease P-E in Pacific/Atlantic as in Broecker (1991)
! [after Oort 1983]: net freshwater loss by Atlantic = 0.32 Sv
! here add/remove total extra1a, extra1b, extra1c Sv of freshwater
! equally by area in Pac/Atl resp.

      do j=1,jmax
         do i=1,imax
            pmeadj(i,j) = 0.  ! P-E adj
         enddo
      enddo

      do j=16,25
         do i=ips(j),ipf(j)        !asurf(j)- grid cell area (m**2)
            pmeadj(i,j) = 1.e6*extra1a/(npac1a*asurf(j))
         enddo
         do i=ias(j),iaf(j)
            pmeadj(i,j) = -1.e6*extra1a/(natl1a*asurf(j))
         enddo
      enddo

      do j=26,51
         do i=ips(j),ipf(j)
            pmeadj(i,j) = 1e6*extra1b/(npac1b*asurf(j))
         enddo
         do i=ias(j),iaf(j)
            pmeadj(i,j) = -1e6*extra1b/(natl1b*asurf(j))
         enddo
      enddo

      do j=52,jmax
         do i=ips(j),ipf(j)
            if(k1(i,j).le.kmax) pmeadj(i,j) = 1e6*extra1c/(npac1c*asurf(j))
         enddo
         do i=ias(j),iaf(j)
            if(k1(i,j).le.kmax) pmeadj(i,j) = -1e6*extra1c/(natl1c*asurf(j))
         enddo
      enddo

! initialize atmosphere

      do j=1,jmax
         do i=1,imax

! initial air temperatures

            tq(1,i,j) = tatm
            tq1(1,i,j) = tq(1,i,j)

! initial specific humidities
! - the ratio of the mass of water vapor in air to the total mass 
! of the mixture of air and water vapor. (grams of vapour per kilogram of air)

! set to relh0_ocean*qsat_ocean over ocean and relh0_land*qsat_atmos over land

            if(k1(i,j).le.kmax)then
               if(ts1(1,i,j,kmax).gt.tsic) then
                  tq(2,i,j) = relh0_ocean*1.e-3*const1* exp(const2*ts1(1,i,j,kmax)/(ts1(1,i,j,kmax)+const3))
               else
                  tq(2,i,j) = relh0_ocean*1.e-3*const1* exp(const4*ts1(1,i,j,kmax)/(ts1(1,i,j,kmax)+const5))
               endif
            else
               if(tq1(1,i,j).gt.0.0) then
                  tq(2,i,j) = relh0_land*1.e-3*const1*exp(const2 *tq1(1,i,j)/(tq1(1,i,j)+const3))
               else
                  tq(2,i,j) = relh0_land*1.e-3*const1*exp(const4 *tq1(1,i,j)/(tq1(1,i,j)+const5))
               endif
            endif

            tq1(2,i,j) = tq(2,i,j)

! other stuff
            qb(i,j) = 0. !heat flux from sea icea into ocean
            evap(i,j) = 0.                         !(should be <0)
            fx0neto(i,j) = 0.

! initialize  sea ice

! thickness ,fractional area and temperature

            varice(1,i,j) = 0.
            varice1(1,i,j) = varice(1,i,j)
            varice(2,i,j) = 0.
            varice1(2,i,j) = varice(2,i,j)
            tice(i,j) = 0.

! rate of change due to thermodynamics

            dtha(1,i,j) = 0.
            dtha(2,i,j) = 0.

! other stuff

            evapsic(i,j) = 0. ! evaporation over sea ice
            fx0sic(i,j) = 0.   !heat flux into sea ice
            fxlw(i,j) = 0.     !net longwave heat flux into atmosphere
		                      ! over open ocean
            fxsen(i,j) = 0. !sensible heat flux from ocean to atmosphere
         enddo
      enddo

! set up runoff catchment data

      call readroff !define runoff matrix for embm version

! diagnostic calculation of global heat source

      ghs = 0.

      end
