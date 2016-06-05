
! invert.f subroutine to invert matrix for barotropic streamfunction
! c grid version variable depth last change
! version for non-trivial coastlines started
! wind part removed to wind.f
! variable drag

      subroutine invert
      include 'var.f90'

      integer i,j,k,l,n,m,ip1,im

      real tv,tv1,rat

      n = imax
      m = jmax + 1
!gap(mpxi*mpxj,2*mpxi+3) - matrix used in streamfunction calculation
      do 710 i=1,n*m
         do 710 j=1,2*n+3
            gap(i,j)=0.
 710  continue

! Set equation at Psi points, assuming periodic b.c. in i.
! Cannot solve at both i=0 and i=imax as periodicity => would
! have singular matrix. At dry points equation is trivial.

      do 650 i=1,imax
         do 650 j=0,jmax
            k=i + j*n   !s - sin of latitude
		  !rh - reciprocal of ocean depth at u,v,T points
            if(max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.kmax)then
               tv = (s(j+1)*rh(1,i,j+1) - s(j)*rh(1,i,j))/ds/dphi2
               tv1 = (sv(j)*rh(2,i+1,j) - sv(j)*rh(2,i,j))/dphi/ds2
          ! c(j) - cos of latitude,  drag(2,maxi+1,maxj) - linear drag
               gap(k,2)= drag(1,i,j)*c(j)*c(j)/ds/ds*rh(1,i,j) + tv1

               l=n+1
! for periodic boundary in i
               if(i.eq.1)l=2*n+1

               gap(k,l) = drag(2,i,j)*rcv(j)*rcv(j)*rdphi*rdphi *rh(2,i,j) - tv

               gap(k,n+2) = - (drag(2,i,j)*rh(2,i,j) + drag(2,i+1,j)*rh(2,i+1,j))&
                            /cv(j)/cv(j)/dphi /dphi - (drag(1,i,j)*c(j)*c(j)*rh(1,i,j)&
                             + drag(1,i,j+1)*c(j+1)*c(j+1)*rh(1,i,j+1))/ds/ds

               l=n+3
! for periodic boundary in i
               if(i.eq.imax)l=3

               gap(k,l)= drag(2,i+1,j)/cv(j)/cv(j)/dphi/dphi*rh(2,i+1,j) + tv
               gap(k,2*n+2)= drag(1,i,j+1)*c(j+1)*c(j+1)/ds/ds*rh(1,i,j+1) - tv1

            else
               gap(k,n+2) = 1
            endif
 650  continue

! now invert the thing

      do 720 i=1,n*m-1
         im=min(i+n+1,n*m)
         do 725 j=i+1,im
            rat=gap(j,n+2-j+i)/gap(i,n+2)
            ratm(j,j-i)=rat
            if(rat.ne.0)then
               do 730 k=n+2-j+i,2*n+3-j+i
                  gap(j,k)=gap(j,k) - rat*gap(i,k+j-i)
 730           continue
            endif
 725     continue
 720  continue

      end