
! wind.f sets wind stress forcing for barotropic streamfunction
! separated from stream.f

      subroutine wind_mean
      include 'var.f90'

      integer i,j,k,l,n,m,ip1

      n = imax
      m = jmax + 1

! calculate constant wind stress part of forcing
! noting that tau not currently defined outside domain
!gb(maxi*maxj)- source term in streamfunction calculation, then solution of
!                             same equation after call of stream
      do 50 i=1,imax
         do 50 j=0,jmax
            ip1=mod(i,imax) + 1
            k = i + j*n
            if(max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.kmax)then
               if(j.eq.jmax.or.j.eq.0)stop 'wind stress not defined outside domain'
!c* wind stress, tau(1,i,j) is the u component at a u point
!c*              tau(2,i,j) is the v component at a v point
               gb(k) = (tau(2,ip1,j)*rh(2,i+1,j) - tau(2,i,j)*rh(2,i,j))&
                       *rdphi*rcv(j) - (tau(1,i,j+1)*c(j+1)*rh(1,i,j+1) - &
                       tau(1,i,j)*c(j)*rh(1,i,j))*rds
            else
               gb(k) = 0.
            endif
            gbold(k) = gb(k)
 50   continue

      end
