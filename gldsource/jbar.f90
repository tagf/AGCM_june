
! jbar.f subroutine to calculate jbar forcing for streamfunction
! variable depth
! adds time dependent pressure torque term to forcing in psi eqn
! altered for non-trivial islands, should give same answers

      subroutine jbar
      include 'var.f90'

      logical getj(maxi,maxj)
      common /lars/getj

      integer i,j,k,l,n,m,km,ip1

      real tv1, tv2, tv3, tv4

      n = imax
      m = jmax + 1

! calculate easy part of double p integral

      do j=1,jmax
         do i=1,imax
! if(mk(i,j).gt.1)then
! bp now defined for k.ge.k1 at all wet points of depth > 1. Other points unset
            if(k1(i,j).lt.kmax)then
               bp(i,j,k1(i,j)) = 0.
               do k=k1(i,j)+1,kmax
                  bp(i,j,k) = bp(i,j,k-1) - (rho(i,j,k) + rho(i,j,k-1)) *dza(k-1)*0.5
               enddo
            endif
         enddo
      enddo

! sbp now defined if mk.gt.0 ie all wet points

      do j=1,jmax
         do i=1,imax
! if(mk(i,j).gt.1)then
            if(mk(i,j).gt.0)then
               sbp(i,j) = 0.
               do k=mk(i,j)+1,kmax
                  sbp(i,j) = sbp(i,j) + bp(i,j,k)*dz(k)
               enddo
            endif
         enddo
      enddo

! periodic b.c. required for Antarctic island integral
! (if bp was defined everywhere would not need ifs)

      do j=1,jmax
         if(k1(1,j).lt.kmax)then
            do k=k1(1,j),kmax
               bp(imax+1,j,k) = bp(1,j,k)
            enddo
         endif
         if(k1(1,j).le.kmax)then
            sbp(imax+1,j) = sbp(1,j)
         endif
      enddo


! calc tricky bits and add to source term for Psi, ip1 copes with
! periodicity in i, j bdy points are ignored assuming no flow out
! of north or south of domain.

      do 70 j=1,jmax-1
         do 70 i=1,imax
            ip1=mod(i,imax) + 1
            l = i + j*n
            if(getj(i,j))then
               tv1 = 0.
               do k=ku(2,ip1,j),mk(ip1,j+1)
                  tv1 = tv1 + bp(ip1,j+1,k)*dz(k)
               enddo
               tv2 = 0.
               do k=ku(2,ip1,j),mk(ip1,j)
                  tv2 = tv2 + bp(ip1,j,k)*dz(k)
               enddo
               tv3 = 0.
               do k=ku(2,i,j),mk(i,j+1)
                  tv3 = tv3 + bp(i,j+1,k)*dz(k)
               enddo
               tv4 = 0.
               do k=ku(2,i,j),mk(i,j)
                  tv4 = tv4 + bp(i,j,k)*dz(k)
               enddo
               gb(l) = gbold(l) + ((tv3 + sbp(i,j+1) - tv4 - sbp(i,j))&
                *rh(2,i,j) - (tv1 + sbp(ip1,j+1) - tv2 - sbp(ip1,j))&
                 *rh(2,ip1,j))*rdphi*rds
               tv1 = 0.
               do k=ku(1,i,j+1),mk(ip1,j+1)
                  tv1 = tv1 + bp(ip1,j+1,k)*dz(k)
               enddo
               tv2 = 0.
               do k=ku(1,i,j),mk(ip1,j)
                  tv2 = tv2 + bp(ip1,j,k)*dz(k)
               enddo
               tv3 = 0.
               do k=ku(1,i,j+1),mk(i,j+1)
                  tv3 = tv3 + bp(i,j+1,k)*dz(k)
               enddo
               tv4 = 0.
               do k=ku(1,i,j),mk(i,j)
                  tv4 = tv4 + bp(i,j,k)*dz(k)
               enddo
               gb(l) = gb(l) + ((tv1 + sbp(ip1,j+1) - tv3 - sbp(i,j+1))&
                *rh(1,i,j+1) - (tv2 + sbp(ip1,j) - tv4 - sbp(i,j))*rh(1,i,j))&
                 *rdphi*rds
            else
               gb(l) = gbold(l)
            endif
   70 continue
! don't currently need the following lines as gb already
! reset to zero at these boundaries by mains
! do i=1,imax
! gb(i) = gbold(i)
! gb(i+jmax*imax) = gbold(i+imax*jmax)
! enddo

      end
