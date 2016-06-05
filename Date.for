      subroutine PrnDate(ndevice)

      INTEGER*2 tmpday, tmpmonth, tmpyear
      integer ndevice
C
C     Show current date .
C
      CALL GETDAT(tmpyear, tmpmonth, tmpday)
      WRITE (ndevice, 1) tmpmonth, tmpday, tmpyear

   1  FORMAT (1X/,' DATE  ',I2, '/', I2.2, '/', I4.4)
      return
      END

