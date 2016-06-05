
! subroutine tstipa.f atmospheric timestep for gldstn
! JGS iterative implicit version
! NRE 2-step version, 2nd attempt
! to recover old explicit code change the lines indicated
! coefficient of implicitness cimp included
! cimp=1 fully implicit, cimp=0 explicit
! coeffs for iterative implicit scheme are defined at cell faces.
! eg flux across east face = cie(i)*T(i+1) + ciw(i)*T(i)
! converted from ocean to atmosphere

      subroutine tstipa

      include 'var.f90'
!	include 'varAGCM.f90' 
	     real  TS_atm_for_oc(0:72+1,0:72+1), QS_atm_for_oc(0:72+1,0:72+1),&
             PREC_atm_for_oc(0:72+1,0:72+1)  
         common /varsAGCM/  TS_atm_for_oc, QS_atm_for_oc, PREC_atm_for_oc

      real tv, ups, ups0, pec, diffpp, cimp, centre, dtloc

      real cie(0:maxi,0:maxj),ciw(0:maxi,0:maxj), cin(0:maxi,0:maxj),cis(0:maxi,0:maxj)
      real tq2(0:maxi+1,0:maxj+1)

! iterations to solve timestep

      integer iits, nii

! implicit
      parameter (nii=4, ups0=999, cimp=0.5)
!      parameter (nii=4, ups0=0.8, cimp=1.0)
! recover old explicit
! parameter (nii=8, ups0=0.0, cimp=0.0)

      integer i, j, l, istep

      logical correct

      parameter(correct=.true. )

      dtloc = dtatm

! set b.c's on local variables

      do i=0,imax
         cin(i,0) = 0.
         cis(i,0) = 0.
         tq2(i,0) = 0.
         cin(i,jmax) = 0.
         cis(i,jmax) = 0.
         tq2(i,jmax+1) = 0.
      enddo

      do l=1,2
         do j=1,jmax
            do i=1,imax
! flux to east
               cie(i,j) = beta(l)*uatm(1,i,j)*rc(j)*0.5*rdphi
               diffpp = diffa(l,1,j) + (2-l)*diffmod0*max(0.0,min(1.0, &
                         (pptn(i,j)-ppmin)/(ppmax-ppmin)))

               tv = rc(j)*rc(j)*rdphi*diffpp*rdphi
! recover old explicit
! ups = sign(ups0, uatm(1,i,j))
               pec = beta(l)*uatm(1,i,j)*dphi/diffpp
               ups = pec / (2.0 + abs(pec))
               ciw(i,j) = cie(i,j)*(1+ups) + tv
               cie(i,j) = cie(i,j)*(1-ups) - tv
! flux to north
               cin(i,j) = cv(j)*beta(l)*uatm(2,i,j)*0.5*rds
               diffpp = diffa(l,2,j) + (2-l)*diffmod0*max(0.0,min(1.0, &
                           (pptn(i,j)-ppmin)/(ppmax-ppmin)))
               tv = cv(j)*cv(j)*rds*diffa(l,2,j)*rds
! recover old explicit
! ups = sign(ups0, uatm(2,i,j))
               pec = beta(l)*uatm(2,i,j)*ds/diffpp
               ups = pec / (2.0 + abs(pec))
               cis(i,j) = cin(i,j)*(1+ups) + tv
               cin(i,j) = cin(i,j)*(1-ups) - tv
            enddo !i
         enddo  !j
         do j=1,jmax
            cie(0,j) = cie(imax,j)
            ciw(0,j) = ciw(imax,j)
         enddo

! iterate to solve timestep

         do iits=1,nii
            do j=1,jmax
               do i=1,imax
                  tq2(i,j) = cimp*tq(l,i,j) + (1.0 - cimp)*tq1(l,i,j)
               enddo
            enddo
            do j=1,jmax
               tq2(0,j) = tq2(imax,j)
               tq2(imax+1,j) = tq2(1,j)
            enddo
            do j=1,jmax
               do i=1,imax
                  centre = dtloc*(ciw(i,j) - cie(i-1,j) + cis(i,j) - cin(i,j-1))
!tqa(2,maxi,maxj) - heat (1) and specific humidity (2) fluxes into atmosphere
                  tq(l,i,j) = (tq1(l,i,j)*(1.0 - (1.0-cimp) *centre) - dtloc*(-tqa(l,i,j)&
                              + cie(i,j)  *tq2(i+1,j) - ciw(i-1,j)*tq2(i-1,j) + cin(i,j) &
                               *tq2(i,j+1) - cis(i,j-1)*tq2(i,j-1)))/ (1. + cimp*centre)
               enddo
            enddo
         enddo !iits
         if(correct)then
            do j=1,jmax
               do i=1,imax
                  tq2(i,j) = 0.5*(tq2(i,j) + cimp*tq(l,i,j) + (1.0 - cimp)*tq1(l,i,j))
               enddo
            enddo
            do j=1,jmax
               tq2(0,j) = tq2(imax,j)
               tq2(imax+1,j) = tq2(1,j)
            enddo
            do j=1,jmax
               do i=1,imax

! explicit and conservative corrector step

                     tq(l,i,j) =  tq1(l,i,j) - dtloc*(-tqa(l,i,j) + cie(i,j) &
                                *tq2(i+1,j) - ciw(i-1,j)*tq2(i-1,j) + cin(i,j)&
                                  *tq2(i,j+1) - cis(i,j-1)*tq2(i,j-1)) - dtloc*tq2(i,j)*&
                                  ( ciw(i,j) - cie(i-1,j) + cis(i,j) - cin(i,j-1) )

           ! write (113,*) i,j,tq(1,i,j)-TS_atm_for_oc(i,j) !GGGGGGGGGG
                    ! tq(1,i,j)=TS_atm_for_oc(i,j) !GGGGGGGGGG
                    ! tq(2,i,j)=QS_atm_for_oc(i,j) !GGGGGGGGGG
! heat divergence in W/m**2 (purely diagnostic)

        transptq(l,i,j) = rhoair*hatmbl(l)*cpa*(tq(l,i,j) - tq1(l,i,j))/(dtatm*rsc/usc)

               enddo !i
            enddo !j
         endif !correct
      enddo  !l

! update tq1

      do j=1,jmax
         do i=1,imax
            do l=1,lmax
! tv = abs(tq1(l,i,j) - tq(l,i,j))
               tq1(l,i,j) = tq(l,i,j)
            enddo
         enddo
      enddo

      end
