
! subroutine to define drag matrix for v2_1 version

      subroutine drgset(adrag,drgf,kmxdrg,jeb)
!jeb - width of equatorial region of maximum drag
      include 'var.f90'

      real tmpdrg(0:maxi,0:maxj)
      real adrag, drgf

      integer i, j, i1, i1p, j1, kloc2, kloc4, kmxdrg, jeb

! calculate drag at psi points (temporary variable)
! first find if there is shallow water (k1>kmxdrg) in
! 2x2 cell neighbourhood then in 4x4 nbhd.
! Increase drag near equator, assuming domain is symmetric

      do j=0,jmax
         do i=0,imax
            kloc2 = max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1))
            kloc4 = k1(i,j)
            do j1=max(0,j-1),min(jmax+1,j+2) !=j-1, j,j+1,j+2
               do i1=i-1,i+2  !=i-1,i,i+1,i+2
                  i1p = 1 + mod(imax + i1-1,imax)
                  kloc4 = max(kloc4,k1(i1p,j1))
               enddo
            enddo
! if(kloc2.gt.kmxdrg.or.j.eq.jmax/2)then
            if(kloc2.gt.kmxdrg.or.abs(j-jmax/2).le.jeb)then
               tmpdrg(i,j) = adrag*drgf*drgf
            else if(kloc4.gt.kmxdrg.or.abs(j-jmax/2).eq.jeb+1)then
               tmpdrg(i,j) = adrag*drgf
            else
               tmpdrg(i,j) = adrag
            endif
         enddo
      enddo

! interpolate to velocity points

      do j=1,jmax
         do i=1,imax
            drag(1,i,j) = 0.5*(tmpdrg(i,j) + tmpdrg(i,j-1))
            drag(2,i,j) = 0.5*(tmpdrg(i,j) + tmpdrg(i-1,j))
         enddo
      enddo

! boundary conditions, assuming no flow out of N or S boundary

      do j=1,jmax
         drag(2,imax+1,j) = drag(2,1,j)
      enddo
! write(101,'(e15.5)')drag

      end