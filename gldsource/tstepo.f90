
! subroutine tstep.f for program gldstn
! flux version fully explicit one step second order variable depth

! version with isoneutral diffusion
! notes; fe/fn/fa in previous versions could have been scalars
! error in isoneutral coeffs corrected
! upstream weighting included

      subroutine tstepo(istep)

      include 'var.f90'

      real tv, ups, ups0, pec

      real fe(maxl), fw(maxl), fn(maxl), fs(maxl,maxi), fa(maxl) , fb(maxl,maxi,maxj), fwsave(maxl)
      integer i, j, k, l, istep

! logical vdt
! parameter(vdt = .false.)
      parameter(ups0 = 0.)

!#ifdef diso
      real tec, scc, tatw, dzrho, tvd, ssmax, rdzrho, slim, tv1
      real dxrho(4), dxts(2,4), dyrho(4), dyts(2,4), dzts(2)
      integer ina, nnp, knp

      parameter(ssmax = 10.)

      scc = ec(2)
      limps = 0
!#endif

! if(mod(istep,100 ).eq.1)
! call velc

      dmax = 0.

! find timestep

! if(vdt)then
! endif

! 2nd order explicit step

! lower boundary fluxes

      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,lmax
               fb(l,i,j) = 0.
  220 continue

      do 100 k=1,kmax

! southern boundary fluxes

         j = 1
         do 230 i=1,imax
            do 230 l=1,lmax
               fs(l,i) = 0.
  230    continue
         do 100 j=1,jmax
! western boundary fluxes
            i = 1
            do 210 l=1,lmax
               if (k.ge.max(k1(imax,j),k1(1,j)))then
! western doorway
! ups = sign(ups0, u(1,imax,j,k))
                  pec = u(1,imax,j,k)*dphi/diff(1)
                  ups = pec / (2.0 + abs(pec))
                  fw(l) = u(1,imax,j,k)*rc(j)*((1.-ups)*ts1(l,1,j,k) + (1.+ups)*ts1(l,imax,j,k))*0.5
                  fw(l) = fw(l) - (ts1(l,1,j,k) - ts1(l,imax,j,k)) *rc2(j)*diff(1)
               else
                  fw(l) = 0.
               endif
               fwsave(l) = fw(l)
  210       continue
            do 100 i=1,imax
               do l=1,lmax
! flux to east
                  if(i.eq.imax)then
! eastern edge(doorway or wall)
                     fe(l) = fwsave(l)
                  elseif(k.lt.max(k1(i,j),k1(i+1,j)))then
                     fe(l) = 0.
                  else
! ups = sign(ups0, u(1,i,j,k))
                     pec = u(1,i,j,k)*dphi/diff(1)
                     ups = pec / (2.0 + abs(pec))
                     fe(l) = u(1,i,j,k)*rc(j)*((1.-ups)*ts1(l,i+1,j,k) + (1.+ups)*ts1(l,i,j,k))*0.5
                     fe(l) = fe(l) - (ts1(l,i+1,j,k) - ts1(l,i,j,k)) *rc2(j)*diff(1)
                  endif
! flux to north
                  if(k.lt.max(k1(i,j),k1(i,j+1)))then
                     fn(l) = 0.
                  else
! ups = sign(ups0, u(2,i,j,k))
                     pec = u(2,i,j,k)*ds/diff(1)
                     ups = pec / (2.0 + abs(pec))
                     fn(l) = cv(j)*u(2,i,j,k)*((1.-ups)*ts1(l,i,j+1,k) + (1.+ups)*ts1(l,i,j,k))*0.5
                     fn(l) = fn(l) - cv2(j)*(ts1(l,i,j+1,k) - ts1(l,i,j,k))*diff(1)
                  endif
! flux above
                  if(k.lt.k1(i,j))then
                     fa(l) = 0.
                  elseif(k.eq.kmax)then
                     fa(l) = ts(l,i,j,kmax+1)
                  else
! ups = sign(ups0, u(3,i,j,k))
                     pec = u(3,i,j,k)*dza(k)/diff(2)
                     ups = pec / (2.0 + abs(pec))
                     fa(l) = u(3,i,j,k)*((1.-ups)*ts1(l,i,j,k+1) + (1.+ups)*ts1(l,i,j,k))*0.5
                     fa(l) = fa(l) - (ts1(l,i,j,k+1) - ts1(l,i,j,k)) *rdza(k)*diff(2)
                  endif
               enddo
!#ifdef diso
! isoneutral diffusion
               if(k.ge.k1(i,j).and.k.lt.kmax)then
                  tatw = 0.5*(ts1(1,i,j,k) + ts1(1,i,j,k+1))
                  tec = - ec(1) - ec(3)*tatw*2 - ec(4)*tatw*tatw*3
                  dzrho = (scc*(ts1(2,i,j,k+1) - ts1(2,i,j,k)) - tec*(ts1(1,i,j,k+1) - ts1(1,i,j,k)))*rdza(k)
                  if(dzrho.lt.-1e-12)then
                     rdzrho = 1.0/dzrho
                     tv1 = 0.0
! tracer loop
                     do knp=0,1
                        do nnp=0,1
                           ina = 1+nnp + 2*knp
                           do l=1,lmax
! phi derivatives
                              dxts(l,ina) = (ts1(l,i+nnp,j,k+knp) - ts1(l,i+nnp-1,j,k+knp))*rc(j)*rdphi
! s-derivatives
                              dyts(l,ina) = (ts1(l,i,j+nnp,k+knp) - ts1(l,i,j+nnp-1,k+knp))*cv(j-1+nnp)*rds
                           enddo
                           dxrho(ina) = scc*dxts(2,ina)-tec*dxts(1,ina)
                           dyrho(ina) = scc*dyts(2,ina)-tec*dyts(1,ina)
! mask each portion of flux separately
                           dxrho(ina) = dxrho(ina)*0.5* (1+sign(1,k+knp-k1(i-1+2*nnp,j)))
                           dyrho(ina) = dyrho(ina)*0.5* (1+sign(1,k+knp-k1(i,j-1+2*nnp)))
! calculate diagonal part
                           tv1 = tv1 + dxrho(ina)*dxrho(ina) + dyrho(ina)*dyrho(ina)
                        enddo
                     enddo
                     tv1 = 0.25*tv1*rdzrho*rdzrho
! limit flux by factor slim for large slope
                     if(tv1.gt.ssmax)then
! slim = ssmax/tv1
                        slim = ssmax*ssmax/(tv1*tv1)
! count flux-limited points
                        limps = limps + 1
                     else
                        slim = 1.0
                     endif
                     tv1 = tv1*slim*diff(1)*rdza(k)
! test vertical diffusion number
                     tv = tv1*dt(k)*rdza(k)
                     if(tv.gt.dmax)then
                        dmax = tv
                     endif
                     do l=1,lmax
                        dzts(l) = (ts1(l,i,j,k+1) - ts1(l,i,j,k))*rdza(k)
! add isoneutral vertical flux
                        tv = 0.
                        do ina=1,4
                           tv = tv + (2.*dzrho*dxts(l,ina) - dxrho(ina)*dzts(l))&
                             *dxrho(ina) + (2.*dzrho*dyts(l,ina) - dyrho(ina)*dzts(l))*dyrho(ina)
                        enddo
                        tv = 0.25*slim*diff(1)*tv/(dzrho*dzrho)
                        fa(l) = fa(l) + tv
! if(abs(tv).gt.1e0)
! 1                     write(6,'(3i4,8e14.4)')i,j,k,fa(l),tv
                     enddo
                  endif
               endif
!#endif
               do l=1,lmax
                  tv = 0.
                  if(k.ge.k1(i,j))then
                     ts(l,i,j,k) = ts1(l,i,j,k) - dt(k)*( - tv + (fe(l) - fw(l))*rdphi&
                              + (fn(l) - fs(l,i))*rds + (fa(l) - fb(l,i,j))*rdz(k))
                  endif

                  fw(l) = fe(l)
                  fs(l,i) = fn(l)
                  fb(l,i,j) = fa(l)
               enddo
               rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k) + ec(3)*ts(1,i,j,k)**2&
                        + ec(4)*ts(1,i,j,k)**3
  100 continue

! convection
      call co(ts)

! periodic b.c. for rho (required at wet points)
! isoneutral code also needs ts1 bc.

      do j=1,jmax
         do k=k1(0,j),kmax
            rho(0,j,k) = rho(imax,j,k)
            do l=1,lmax
               ts1(l,0,j,k) = ts(l,imax,j,k)
            enddo
         enddo
         do k=k1(imax+1,j),kmax
            rho(imax+1,j,k) = rho(1,j,k)
            do l=1,lmax
               ts1(l,imax+1,j,k) = ts(l,1,j,k)
            enddo
         enddo
      enddo

      do 10 k=1,kmax
         do 10 j=1,jmax
            do 10 i=1,imax
               do 10 l=1,lmax
                  if(k.ge.k1(i,j)) ts1(l,i,j,k) = ts(l,i,j,k)
   10 continue

! t = istep*dt(kmax) + t0

      end
