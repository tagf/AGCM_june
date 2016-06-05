
! subroutine surflux.f for program gldstn
! computes heat and freshwater fluxes at atmosphere-ocean/seaice/land interface

      subroutine surflux

      include 'var.f90'

      real ce, ch, cesic, chsic, rq, tv0,tv1,tv2,tv3, tol, set0
      real albsic, fxswsic , fxlwsic, fxsensic , fx0oa, fx0sica
      real qsatsic, zeroc, alw, ticold, cfxsensic, salt, dho, dhsic
      real dhadj,cetemp  !cetemp - for pollution
!##### debug
      real amin,amax,sum1,sum2,sum3,vsc
      integer iamin,iamax,jamin,jamax
 !##########
      parameter(zeroc = 273.15)

      integer i, j, l, iter, itice,i1, j1
      integer YearN1, YearN !for CO2 up to year 2100

!c	integer oilx1, oilx2, oily1, oily2

      parameter(itice = 20 , tol=1.e-10  )

! path of oil
!c	oilx1 = 17
!c	oilx2 = 20
!c	oily1 = 25
!c	oily2 = 27


! initialize integrated runoff
      do i=1,imax
         do j=1,jmax
            runoff(i,j) = 0.
         enddo
      enddo
      do i=1,imax
         do j=1,jmax
!c****           write (555,'(2i3,e12.3)') i,j, tq1(1,i,j)
         enddo
      enddo

!****** for global heat source diagnostic (sb diaga.f)
      ghs=0.
!####### CO2 up to 2100
!      YearN=istepCO2/365+1980
!      YearN1=istepCO2/3650+1 !10 year counter
!     if (YearN1.lt.13) then
!       CO2=CO2+(A2_CO2(YearN1+1)-A2_CO2(YearN1))/3650.
       if (mod(istepCO2,365)==1) then
        write (6,*) ' CO2: ',istepCO2,YearN, CO2(1,1)
       endif
!      endif

! main i,j loop to compute surface flux terms
!c***      print*,'sur',const1,const5
      do i=1,imax
         do j=1,jmax

! pptn (over land and ocean) - precipitation
! - need saturation vapour pressure, relative humidity
!c***      print*,'sur',i,j,tq1(1,i,j),qsata(i,j)
!   tq - temp (1) and specific humidity (2) in atmosphere
!   tq1 - same at previous timestep
!qsata - saturation specific humidity in atm. for precipitation
            qsata(i,j) = 1.e-3*const1*exp(const4*tq1(1,i,j)/(tq1(1,i,j)+const5))

! relative humidity
! - the ratio of the mass of water vapor in air to the total mass 
! of the mixture of air and water vapor. (grams of vapour per kilogram of air)
! nre remove this division??
            rq = tq1(2,i,j)/qsata(i,j)
! rmax=0.85 - threshold relative humidity above which precipitation(m/s) occurs
! rhoao - ratio of air to ocean density
!hatmbl(2) - atmospheric boundary layer height (m) for heat (1) and humidity (2)
!dt(maxk) -ocean and sea-ice timestep, variable in z option never used
           pptn(i,j) = max(0.0,(tq1(2,i,j) - rmax*qsata(i,j))*rhoao*hatmbl(2)/(tsc*dt(kmax)))

! nre instantaneous precipitation

            tq1(2,i,j) = min(tq1(2,i,j),rmax*qsata(i,j))

! use climatological albedo in calculating incoming shortwave radiation

            albedo(i,j) = albcl(i,j)
            !##### snow albedo
       !     if (tq1(1,i,j).lt.-5..and.k1(i,j).gt.kmax) then
       !       albedo(i,j)= 0.85
       !     endif

! shortwave radiation  - flux of shortwave radiation:

            fxsw(i,j) = solfor(j)*(1. - albedo(i,j))

! outgoing planetary longwave

            tv0 = b00 + rq*(b10 + b20*rq)
            tv1 = b01 + rq*(b11 + b21*rq)
            tv2 = b02 + rq*(b12 + b22*rq)
            tv3 = b03 + rq*(b13 + b23*rq)

! update co2 concentration should be moved to mains or tstepa??
! compound increase at fractional rate per dtatm, rate_co2 (see gseta.f)

   !old         co2(i,j) = co2(i,j) + rate_co2*co2(i,j)
		   if (co2(i,j) + delta_co2<coef_co2*co20) then
				co2(i,j) = co2(i,j) + delta_co2
	     end if
!fxplw(maxi,maxj) -flux of ('planetary') longwave radiation to outer space
            fxplw(i,j) = tv0 + tq1(1,i,j)*(tv1 + tq1(1,i,j)*(tv2 + tq1(1,i,j)*tv3))- delf2x*log(co2(i,j)/co20)
! latent heat flux into atmos associated with condensation
! nre with no account taken of snow melting, must assume all pptn is rain
! hlv = 2.501e6  - latent heat of vapourisation (J/kg)
            fxlata(i,j) = rho0*pptn(i,j)*hlv

!c-----------------------------------------------------------------------
! calculate terms over ocean or ice
!c-----------------------------------------------------------------------
            if(k1(i,j).le.kmax) then

! longwave radiation

               alw = tq1(1,i,j)+zeroc    !zeroc = 273.15
               alw = alw * alw
               alw = alw * alw
               alw = ema * alw  !ema = 0.85 * sigma;  sigma = 5.67e-8

! surface salinity-dependent freezing point:

               salt = saln0+ts1(2,i,j,kmax)
! salt = ts1(2,i,j,kmax)
!             freeze temperature
               tsfreez(i,j) = salt*(-0.0575 + 0.0017*sqrt(salt)- 0.0002*salt)
! or constant:  tsic = -1.8
! tsfreez(i,j) = tsic

! maximum amount of heat available in first layer
! nre rsictscsf must be changed if dt>17.5 days, see gseta
!     qb(i,j) - heat flux from sea icea into ocean (should be <0)
               qb(i,j) = rsictscsf*(tsfreez(i,j)-ts1(1,i,j,kmax))


!c-----------------------------------------------------------------------
! calculate terms over ice
!c-----------------------------------------------------------------------
               if(varice1(2,i,j).gt.0.0)then

! let albedo over sea ice vary as a function of tair (Holland et al. 1993)

                  albsic = max(0.20,min(0.6,0.40 - 0.04*tq1(1,i,j)))
                  albedo(i,j)= albsic !#####
                  fxswsic = solfor(j)*(1. - albsic)

! first need to calculate T_ice - iterations
                  do iter=1,itice   !=20
                     ticold = tice(i,j)
! Dalton number
                     cesic = 1.0e-3*(1.0022 - 0.0822*(tq1(1,i,j)- ticold) + 0.0266*usurf(i,j))
                     cesic = max(6.0e-5,min(2.19e-3,cesic))

                     chsic = 0.94*cesic

! sensible heat flux
                     cfxsensic = rhoair*chsic*cpa*usurf(i,j)

                     qsatsic = 1e-3*const1*exp(const2*ticold/(ticold + const3))

                     evapsic(i,j) = (qsatsic - tq1(2,i,j))*rhoao*cesic*usurf(i,j)
!tice(maxi,maxj)   - surface temperature of sea ice
!varice1 - sea ice variables: average height (1) and fractional area (2)
                     tice(i,j) = (difsic*tsfreez(i,j) + varice1(1,i,j)*&
                     (cfxsensic*tq1(1,i,j)+ (1-ca(i,j))*fxswsic + alw- &
                     rho0*hls*evapsic(i,j) + emo*(ticold+zeroc)**3*&
                     (3.0*ticold-zeroc)))/(difsic + varice1(1,i,j)*&
                     (cfxsensic + 4.0*emo*(ticold+zeroc)**3))
!hls = 2.835e6 - latent heat of sublimation
	               !### my change - for convergence
                     if(tice(i,j).gt.tfreez) tice(i,j)=tfreez
                     if(abs(tice(i,j) - ticold).lt.tol) goto 10
                  enddo  !iter
! print*,'tice non-convergence at',
! 1                    i,j,tice(i,j),tice(i,j)-ticold

   10             tice(i,j) = min(tfreez,tice(i,j))

                  fxlwsic = emo*(tice(i,j)+zeroc )**4 - alw

                  fxsensic = cfxsensic*(tice(i,j) - tq1(1,i,j))

                  fx0sic(i,j) = (1-ca(i,j))*fxswsic -fxsensic- fxlwsic - rho0*hls*evapsic(i,j)

                  fx0sica = ca(i,j)*fxswsic + fxlata(i,j)+ fxsensic + fxlwsic- fxplw(i,j)
                           !rrholf - reciprocal of rho*Lf
                  dhsic = rrholf*(qb(i,j) - fx0sic(i,j)) - rhooi*evapsic(i,j)
               else
                  fx0sica = 0.0
                  dhsic = 0.
                  evapsic(i,j) = 0.
                  tice(i,j) = 0.
               endif

!c-----------------------------------------------------------------------
! over open ocean
!c-----------------------------------------------------------------------

               fxlw(i,j) = emo*(ts1(1,i,j,kmax)+zeroc)**4 - alw

! Dalton number

               ce = 1.0e-3*(1.0022 - 0.0822*(tq1(1,i,j)-ts1(1,i,j,kmax)) + 0.0266*usurf(i,j))
               ce = max(6.0e-5,min(2.19e-3,ce))
	! ######### Pacific ocean pollution exper. 2*evapor
!	         cetemp=ce
!	         if (j==26.and.(i==4.or.i==5.or.i==6.or.i==12.or.i==13))
!     *			  ce=cetemp*2.
!	         if (j==27.and.(i==4.or.i==5.or.i==6.or.i==7.or.i==10
!     *			 .or.i==11.or.i==12.or.i==13.or.i==14)) ce=cetemp*2.
!	         if (j==28.and.(i==5.or.i==6.or.i==7.or.i==8.or.i==10
!     *			 .or.i==11.or.i==12.or.i==13.or.i==14)) ce=cetemp*2.
!	         if (j==29.and.(i==5.or.i==6.or.i==7.or.i==8.or.i==10
!     *			 .or.i==11.or.i==12.or.i==13)) ce=cetemp*2.
!	         if (j==30.and.(i==6.or.i==7)) ce=cetemp*2.
	! #########
               ch = 0.94*ce

! sensible heat flux from ocean to atmosphere
! cpa = 1004.  !specific heat capacity of air (J/kg/Ê)
	       fxsen(i,j) = rhoair*ch*cpa*usurf(i,j)*(ts1(1,i,j,kmax)-tq1(1,i,j))

! evaporation/sublimation rate
!qsato(maxi,maxj)  - saturation specific humidity over ocean for evaporation
               qsato(i,j) = 1.e-3*const1*exp(const4*ts1(1,i,j,kmax)/(ts1(1,i,j,kmax)+const5))

               evap(i,j) = (qsato(i,j) - tq1(2,i,j))*rhoao*ce*usurf(i,j)

! decrease of evaporation in spot of oil		
!c			open (56,file='C:\c-goldstein\results\results_oil')
!c			if 	((i.ge.oilx1).and.(i.le.oilx2).and.(j.ge.oily1)
!c	1			.and.(j.le.oily2)) 	then
!c				write(56, '(2i2,e15.5)') i, j, evap(i,j)
!c				evap(i,j) = evap(i,j)*0.5
!c				evap(i,j) = 0
!c				write(56, '(2i2,e12.3)') i, j, evap(i,j)
!c			end if


! net heat flux into atmosphere

               fx0oa = ca(i,j)*fxsw(i,j) + fxlata(i,j) + fxsen(i,j) + fxlw(i,j) - fxplw(i,j)

! add proportions over open ocean and sea ice

               fx0a(i,j) = (1.-varice1(2,i,j))*fx0oa+ varice1(2,i,j)*fx0sica

! heat flux from atmosphere into open ocean
! hlv = 2.501e6  - latent heat of vapourisation (J/kg)
               fx0o(i,j) = (1.-ca(i,j))*fxsw(i,j) - fxsen(i,j)- fxlw(i,j) - rho0*hlv*evap(i,j)

! net heat flux into ocean from atmosphere and sea ice
! including possible ice growth over open ocean
!varice1 - sea ice variables: average height (1) and fractional area (2)

               fx0neto(i,j) = varice1(2,i,j)*qb(i,j)+ (1.-varice1(2,i,j))*max(qb(i,j),fx0o(i,j))

               dho = max(0.0,rrholf*(qb(i,j) - fx0o(i,j)))

        !  rate of change of sea-ice height and area respectively
               dtha(1,i,j) = varice1(2,i,j)*dhsic + (1-varice1(2,i,j))*dho
               dtha(2,i,j) = max(0.0,rhmin*dho*(1.-varice1(2,i,j)))
               if(varice1(1,i,j).gt.1.e-12) dtha(2,i,j) = dtha(2,i,j)+&
                   min(0.0,0.5*varice1(2,i,j)*dtha(1,i,j)/varice1(1,i,j))

! global heat source diagnostic

               ghs = ghs + varice1(2,i,j)*(fx0sic(i,j) + fx0sica)+ &
               (1.-varice1(2,i,j))*(fx0o(i,j) + fx0oa)+ rhoice*hlf*&
               (dtha(1,i,j)+ evapsic(i,j)*varice1(2,i,j)*rhooi)
            else
!c-----------------------------------------------------------------------
! calculate terms over land
!c-----------------------------------------------------------------------
! 'zero's' set in gseta. Evap=0 over land (no moisture or heat capacity).
       !fx0a - net heat flux into atmosphere
               fx0a(i,j) = fxsw(i,j) + fxlata(i,j)- fxplw(i,j)

! runoff scheme: in the case of land pptn .ne. zero, find nearest ocean
! gridbox/gridboxes and add as runoff there

               if(pptn(i,j).ne.0.0) then
                  runoff(iroff(i,j),jroff(i,j)) = runoff(iroff(i,j),jroff(i,j)) + pptn(i,j)
               endif

               ghs = ghs + fx0a(i,j)
            endif

! perturb runoff by extra0 Sv
! nre extra0 def'n altered to avoid unnecessary division

! ...for Stefan's intercomparison expts (NB. no compenstion elsewhere)...
! over 62 gridboxes south of convective regions (20-50N) in N.Atlantic:
!	     if(j.ge.25.and.j.le.32.and.i.ge.ias(j).and.i.le.iaf(j))
!     1         runoff(i,j) = runoff(i,j) + extra0*rextra0
! or over 18 gridboxes within convective regions (50-70N) in N.Atlantic:
            if(j.ge.66.and.j.le.70.and.i.ge.ias(j).and.i.le.iaf(j))&
            runoff(i,j) = runoff(i,j) + extra0*rextra0

! end of main i,j loop

         enddo !j
      enddo  !i
	if (istepT.le.-2000)then
	open(229,file='test_ts1')
	do i=1,imax
	   do j=1,jmax
	   write (229,101) i,j,ts1(2,i,j,kmax)
         enddo
      enddo
101	format (i3,',',i3,',',e10.2)
	endif
! short i,j loop to add in effects of runoff (non-local in i,j)
! and set up fluxes

      do i=1,imax
         do j=1,jmax

! freshwater forcing in m/s open ocean P-E over ocean gridboxes
! evap zero over land,
! ??nre P goes straight through as no snow allowed, but E is E total
! for atm model and soley E_ocean for the fwflux

            pme(i,j) = pptn(i,j) - evap(i,j)*(1-varice1(2,i,j))- evapsic(i,j)*varice1(2,i,j)

! non-dimensionalize surface fluxes for use in tstepa:

            tqa(1,i,j) = fx0a(i,j)*rfluxsca
! tqa(2,i,j) = - pme(i,j)*rpmesca
! nre instantaneous precipitation
            tqa(2,i,j) = (evap(i,j)*(1-varice1(2,i,j))+ evapsic(i,j)*varice1(2,i,j))*rpmesca

! optional split timestep for atmosphere
! tq(1,i,j) = tq1(1,i,j) + dt(kmax)*tqa(1,i,j)
! tq(2,i,j) = tq1(2,i,j) + dt(kmax)*tqa(2,i,j)
! tq1(1,i,j) = tq(1,i,j)
! tq1(2,i,j) = tq(2,i,j)
! reset tqa to zero if used here (tqa subsequently used in tstepa)
! tqa(1,i,j) = 0.
! tqa(2,i,j) = 0.

! latent heat flux associated with evaporation or sublimation
! purely diagnostic variable

            fxlato(i,j) = rho0*evap(i,j)*hlv*(1-varice1(2,i,j))+ &
                       rho0*hls*evapsic(i,j)*varice1(2,i,j)
         enddo
      enddo

! update ice

      call tstepsic

! modify heat and salinity fluxes into ocean according to sea-ice update;
! prevent <0%, >100% fractional area A and set A=H=0 if H<Hmin
! in which case add an amount -H of ice, hence need to add appropriate
! heat and freshwater fluxes.

      do j=1,jmax
         do i=1,imax
            if(kmax.ge.k1(i,j))then
!         fwfxneto - net freshwater flux into ocean (m/s)
!pme(maxi,maxj)    ; precipitation (P) - evaporation (E) (m/s)
!dtha(2,maxi,maxj) ; rate of change of sea-ice height (1) and area respectively
!rhoio,rhooi - ratio of sea-ice to ocean and ocean to sea-ice densities resp.
! hlf = 3.34e5 - latent heat of fusion of ice (J/kg)
                  !pmeadj(i,j)- P-E adj - implicit fresh water flux
               fwfxneto(i,j) = pme(i,j) + pmeadj(i,j)+ runoff(i,j) - rhoio*dtha(1,i,j)
	!varice - sea ice variables: average height (1) and fractional area (2)
               varice(2,i,j) = max(0.0,min(1.0,varice(2,i,j)))
               if(varice(1,i,j).lt.hmin)then
                  fx0neto(i,j) = fx0neto(i,j) - varice(1,i,j)*rhoice*hlf/(tsc*dt(kmax))
                  fwfxneto(i,j) = fwfxneto(i,j) + varice(1,i,j)*rhoio /(tsc*dt(kmax))
                  ghs = ghs - varice(1,i,j)*rhoice*hlf/(tsc*dt(kmax))
                  do l=1,2
                     varice(l,i,j) = 0.
                  enddo
               endif

! upper boundary conditions for ocean: k=kmax+1
! for v3_1 (implicit) ocean code the source is stored in ts(...,kmax+1)
! non-dimensionalize surface fluxes for use in tstepo
!rfluxsc,rfluxsca - reciprocal dimensional scale values for heat fluxes into ocean and atm.
!rpmesco = rsc*saln0/(dsc*usc) - reciprocal scaling for freshwater forcing of ocean
               ts(1,i,j,kmax+1) = - fx0neto(i,j)*rfluxsc
               ts(2,i,j,kmax+1) = fwfxneto(i,j)*rpmesco
               ts1(1,i,j,kmax+1) = ts(1,i,j,kmax+1)
               ts1(2,i,j,kmax+1) = ts(2,i,j,kmax+1)

               do l=1,2
                  varice1(l,i,j) = varice(l,i,j)
               enddo
            endif
         enddo
      enddo

      end
