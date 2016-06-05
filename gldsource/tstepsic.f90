
! subroutine tstepsic.f for program gldstn
! updates sea-ice height and area

      subroutine tstepsic

      include 'var.f90'

      integer i, j, l
      real tv

      real fe(2), fw(2), fn(2), fs(2,maxi) ,fwsave(2)

! old code for split timestep
! do i=1,imax
! do j=1,jmax
! if(k1(i,j).le.kmax) then
! do l=1,2
! varice(l,i,j) = varice1(l,i,j)
! 1                          + tsc*dt(kmax)*dtha(l,i,j)
! enddo
! endif
! enddo
! enddo

! 2nd order explicit transport code using upper level ocean velocities

! southern boundary fluxes

      j = 1
      do 230 i=1,imax
         do 230 l=1,2
            fs(l,i) = 0
  230    continue
      do 100 j=1,jmax
! western boundary fluxes
         i = 1
         do 210 l=1,2
            if (kmax.ge.max(k1(imax,j),k1(1,j)))then
! western doorway
               fw(l) = u(1,imax,j,kmax)*rc(j)*(varice1(l,1,j) + varice1(l,imax,j))*0.5
               fw(l) = fw(l) - (varice1(l,1,j) - varice1(l,imax,j)) *rc(j)*rc(j)*rdphi*diff(1)
            else
               fw(l) = 0
            endif
            fwsave(l) = fw(l)
  210       continue
         do 100 i=1,imax
            do 120 l=1,2
! flux to east
               if(i.eq.imax)then
! eastern edge(doorway or wall)
                  fe(l) = fwsave(l)
               elseif(kmax.lt.max(k1(i,j),k1(i+1,j)))then
                  fe(l) = 0.
               else
                  fe(l) = u(1,i,j,kmax)*rc(j)*(varice1(l,i+1,j) + varice1(l,i,j))*0.5
                  fe(l) = fe(l) - (varice1(l,i+1,j) - varice1(l,i,j)) *rc(j)*rc(j)*rdphi*diff(1)
               endif
! flux to north
               if(kmax.lt.max(k1(i,j),k1(i,j+1)))then
                  fn(l) = 0.
               else
                  fn(l) = cv(j)*u(2,i,j,kmax)*(varice1(l,i,j+1) + varice1(l,i,j))*0.5
                  fn(l) = fn(l) - cv(j)*cv(j)*(varice1(l,i,j+1) - varice1(l,i,j))*rds*diff(1)
               endif

               if(kmax.ge.k1(i,j))then
                  varice(l,i,j) = varice1(l,i,j) - dt(kmax)*( (fe(l) - fw(l))*rdphi + (fn(l) - fs(l,i))*rds) + tsc*dt(kmax)*dtha(l,i,j)
               endif
               fw(l) = fe(l)
               fs(l,i) = fn(l)
  120       continue
  100 continue

      end
