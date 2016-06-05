
! ubarsolv.f subroutine to calculate barotropic velocity on c grid
! variable depth last change
! adds time dependent pressure torque term to forcing in psi eqn
! domain extended

      subroutine ubarsolv(ubloc,psiloc)
      include 'var.f90'

      integer i,j,k,l,n,m,km,im

      real ubloc(2,0:maxi+1,0:maxj) !barotropic velocity components
      real psiloc(0:maxi,0:maxj)  !barotropic streamfunction

      n = imax
      m = jmax + 1

! solve Psi equation

      do 10 i=1,n*m-1
         im=min(i+n+1,n*m)
         do 10 j=i+1,im
! ratm(maxi*maxj,maxi+1)- matrix used to invert for the barotropic streamfunction
!gb(maxi*maxj)- source term in streamfunction calculation, then solution of
!                             same equation after call of stream
            gb(j)=gb(j) - ratm(j,j-i)*gb(i)
   10 continue
 !gap(mpxi*mpxj,2*mpxi+3) - matrix used in streamfunction calculation
      gb(n*m)=gb(n*m)/gap(n*m,n+2)
      do 20 i=n*m-1,1,-1
         km=min(n+1,n*m-i)
         do 30 k=1,km
            gb(i)=gb(i) - gap(i,n+2+k)*gb(i+k)
   30    continue
         gb(i)=gb(i)/gap(i,n+2)
   20 continue

! write to Psi for convenience (useful if ACC)

      do j=0,jmax
         do i=1,imax
            k = i + j*n
            psiloc(i,j) = gb(k)
         enddo
         psiloc(0,j) = psiloc(imax,j)
      enddo

! calculate barotropic vely where Psi (and ub) defined

      do j=1,jmax
         do i=1,imax
            ubloc(1,i,j) = -rh(1,i,j)*c(j)*(psiloc(i,j) - psiloc(i,j-1)) *rds
         enddo
      enddo

      do j=1,jmax-1
         do i=1,imax
            ubloc(2,i,j) = rh(2,i,j)*(psiloc(i,j) - psiloc(i-1,j)) *rcv(j)*rdphi
         enddo
      enddo

! set velocity to zero at N and S boundaries
      do i=1,imax
         ubloc(2,i,jmax) = 0.
         ubloc(2,i,0) = 0.
      enddo

! periodic b.c. for ub(2) required only for island integral

      do j=1,jmax
         ubloc(2,imax+1,j) = ubloc(2,1,j)
         ubloc(1,0,j) = ubloc(1,imax,j)
         ubloc(1,imax+1,j) = ubloc(1,1,j)
         ubloc(2,0,j) = ubloc(2,imax,j)
      enddo
      ubloc(2,imax+1,0) = ubloc(2,1,0)
      ubloc(2,0,0) = ubloc(2,imax,0)

      end
