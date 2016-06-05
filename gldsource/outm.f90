! subroutine outm.f writes out data for gldstn
! expanded to write out atmos and sea ice data
! write out: ts(l,i,j,k), u(l,i,j,k),tq(l,i,j),varice(l,i,j),tice(i,j),t
      subroutine outm(unit)

      include 'var.f90'

      integer i, j, k, l, unit

      do 20 j=1,jmax
         do 20 i=1,imax
            do 20 k=1,kmax
               do 30 l=1,lmax
                  if(k.ge.k1(i,j))then
                     write(unit,* )ts(l,i,j,k)
                  else
                     write(unit,* )0.0
                  endif
   30          continue
               do l=1,2
                  write(unit,* )u(l,i,j,k)
               enddo
   20 continue

! EMBM

      do 120 j=1,jmax
         do 120 i=1,imax
            do 120 l=1,2
               write(unit,* )tq(l,i,j)
  120 continue
      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,2
               write(unit,* )varice(l,i,j)
  220 continue

! EMBM for exact continuation need

      do j=1,jmax
         do i=1,imax
            write(unit,* )tice(i,j)
         enddo
      enddo
                    !##############
      write(unit,*)t,sdedy,sdeyr
   10 format(10f10.4)
      end
