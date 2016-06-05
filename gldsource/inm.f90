! subroutine inm.f reads in data for gldstn
! expanded to read in atmos and sea ice data

      subroutine inm(unit)

      include 'var.f90'

      integer i, j, k, l, unit
! read (unit,*)imax,jmax,kmax
!ocean temp (1) (grad C) and salinity (2) array
      read (unit,*)((((ts(l,i,j,k),l=1,lmax),(u1(l,i,j,k),l=1,2), k=1,kmax),i=1,imax),j=1,jmax)
         do j=1,jmax
           do i=1,imax  
            write (1444,*) i,j,ts(1,i,j,kmax),k1(i,j) !GGGGGGGGGG
           enddo
         enddo
   write(1445,'(74F10.4)') ((ts(1,i,j,kmax),i=0,73),j=73,0,-1)!GGGGGGGGGG
 !      stop 'ts  inm'  !test 

! extra read statement for embm atmos
! temp (1) and specific humidity (2) in atmosphere
      read (unit,*)(((tq(l,i,j),l=1,2),i=1,imax),j=1,jmax)

! extra read statement for sea ice
!sea ice variables: average height (1) and fractional area (2)
      read (unit,*)(((varice(l,i,j),l=1,2),i=1,imax),j=1,jmax)

! extra read statement for exact continuation
!  surface temperature of sea ice (grad C)
      read (unit,*)((tice(i,j),i=1,imax),j=1,jmax)
                               !#########
      read (unit,*,end=10) t0, sdedy, sdeyr
      t = t0  !initial time
      print*,'sb inm: t = ',t,' sdedy = ', sdedy,' sdeyr = ',sdeyr
 10   continue
 
      return
      end
