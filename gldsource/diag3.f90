
! diag3.f frequent atmos diagnostics - mean glob air temp, hum,
!                      hice, hice*A, total number of ice points

      subroutine diag3(sum1,sum2,sum3,sum4,iice)

      include 'var.f90'

      real sum1, sum2, sum3, sum4, vol1, vol2

      integer i,j,iice

      sum1=0.
      sum2=0.
      vol1=0.
      vol2=0.

      do j=1,jmax
         do i=1,imax
            sum1 = sum1 + tq_avr(1,i,j)*dphi*ds*hatmbl(1)
	    vol1 = vol1 + dphi*ds*hatmbl(1)
            sum2 = sum2 + tq_avr(2,i,j)*dphi*ds*hatmbl(2)
	    vol2 = vol2 + dphi*ds*hatmbl(2)
         enddo
      enddo

      sum1 = sum1/vol1 !temp
      sum2 = sum2/vol2 !hum
	!########################
      sum3 = 0.
      sum4 = 0.
	iice = 0
      do j=1,jmax
         do i=1,imax
	         if (ice_avr(2,i,j).GT.0.) then
			  sum3 = sum3 +ice_avr(1,i,j) !hice
			  sum4 = sum4 +ice_avr(1,i,j)/ice_avr(2,i,j)
	          iice = iice+1  !ice points
               endif
         enddo
      enddo
	if (iice==0) iice=1
	sum3 = sum3/iice
	sum4 = sum4/iice
      !#########################
      end
