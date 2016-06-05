
! subroutine tstepa.f for program gldstn
! transports tair, qair meridionally and vertically
! updates tair, qair in lower atmosphere

! flux version fully explicit one step second order variable depth

      subroutine tstepa

      include 'var.f90'

      real tv

      real fe(2), fw(2), fn(2), fs(2,maxi), fa(2), fb(2,maxi,maxj) ,fwsave(2)

      real diffextra

      integer i, j, l

! 2nd order explicit step

! lower boundary fluxes

      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,2
               fb(l,i,j) = 0.
  220 continue

! southern boundary fluxes

         j = 1
         do 230 i=1,imax
            do 230 l=1,2
               fs(l,i) = 0.
  230    continue

         do 100 j=1,jmax
! western boundary fluxes
            i = 1
            do 210 l=1,2
! western doorway
!uatm(2,maxi,maxj) - prescribed atmospheric advective velocities
               fw(l) = beta(l)*uatm(1,imax,j)*rc(j)*(tq1(l,1,j) + tq1(l,imax,j))*0.5
! add zonal heat diffusion
               diffextra = (2-l)*diffmod0*max(0.0,min(1.0, (pptn(i,j)-ppmin)/(ppmax-ppmin)))
!diffa(2,2,maxj) - atmospheric diffusivity (l,m,j) where l denotes T or Q, m denotes
!                  u or v points, for each j point; rc = 1/cos(j)
               fw(l) = fw(l) - (tq1(l,1,j) - tq1(l,imax,j)) *rc(j)*rc(j)*rdphi*(diffa(l,1,j)+diffextra)
               fwsave(l) = fw(l)
210         continue

            do 100 i=1,imax
               do 120 l=1,2
! flux to east
                  if(i.eq.imax)then
! eastern edge(doorway or wall)
                     fe(l) = fwsave(l)
                  else
                     fe(l) = beta(l)*uatm(1,i,j)*rc(j)*(tq1(l,i+1,j) + tq1(l,i,j))*0.5
! add zonal heat diffusion
                     diffextra = (2-l)*diffmod0*max(0.0,min(1.0, (pptn(i,j)-ppmin)/(ppmax-ppmin)))
                     fe(l) = fe(l) - (tq1(l,i+1,j) - tq1(l,i,j)) *rc(j)*rc(j)*rdphi*(diffa(l,1,j)+diffextra)
                  endif
! flux to north
                  if(j.ne.jmax)then
! except northermost gridpoints

                     fn(l) = cv(j)*beta(l)*uatm(2,i,j)*(tq1(l,i,j+1) + tq1(l,i,j))*0.5
! add meridional heat diffusion
                     diffextra = (2-l)*diffmod0*max(0.0,min(1.0, (pptn(i,j)-ppmin)/(ppmax-ppmin)))
                     fn(l) = fn(l) - cv(j)*cv(j)*(diffa(l,2,j) + diffextra) *(tq1(l,i,j+1) - tq1(l,i,j))*rds
                  else
                     fn(l) = 0.
                  endif

! if(m.eq.2)then
! zero flux through top of upper layer
! fa(l) = 0
! else
! surface flux
!tqa(2,maxi,maxj) - heat (1) and specific humidity (2) fluxes into atmosphere
                     fb(l,i,j) = tqa(l,i,j)
! flux through interface between lower and upper atmosphere
                     fa(l) = 0.
! endif
!c****	print *,i,j,tq(l,i,j),fe(l),fw(l),fn(l),fs(l,i),fa(l),fb(l,i,j)
      if (l==1)  then
	  write (555,'(2i3,7e12.3)') i,j,tq(l,i,j),fe(l),fw(l),fn(l), &
	                                fs(l,i),fa(l),fb(l,i,j)
	 endif
                     tq(l,i,j) = tq1(l,i,j) - dtatm*( (fe(l) - fw(l))*rdphi + &
                                 (fn(l) - fs(l,i))*rds + fa(l) - fb(l,i,j))
! 3                           + (fa(l) - fb(l,i,j))/(hatmbl(l)/dsc))

! nre dimless height of atm set to 1, implies factor of h in fluxsc

! heat divergence in W/m**2 (purely diagnostic)

        transptq(l,i,j) = rhoair*hatmbl(l)*cpa*(tq(l,i,j) - tq1(l,i,j))/(dtatm*rsc/usc)

                  fw(l) = fe(l)
                  fs(l,i) = fn(l)
                  fb(l,i,j) = fa(l)
  120          continue

  100 continue

      do 10 j=1,jmax
         do 10 i=1,imax
            do 10 l=1,2
               tq1(l,i,j) = tq(l,i,j)
   10 continue

      end
