subroutine smooth(PROP)

!  This subroutine smooths the variable "PROP" (i.e. it adds the
!  artificial viscosity) by taking (1-SF) x the calculated value of
!  "PROP" + SF x (The average of the surrounding values of "PROP").
!  Where SF is the smoothing factor.

	use common_block
!
       implicit none
!
       integer ( kind = 4 ) :: i,j,k,IP1,IM1
       real (kind = dp ) :: SF,SFM1,aux,aux2,aux3,aux4
       real (kind = dp ) :: AVG,AVG1,AVGNJ,AVG2,AVGNK
       real (kind = dp ) :: PROP(NI,NJ,NK),STORE(NI,NJ,NK)
!
!   To avoid using already smoothed values for smoothing other values
!   the smoothed values are initially stored in an array "STORE".
!
       SF   = SMOOTH_FAC
      SFM1 = 1.0 - SF
      aux = 1/6.
      aux2 = 1/5.
      aux3 = 1/4.
      aux4 = 1/3.
!
       do i=1,NI
	   IP1  = i+1
	   IF(i.EQ.NI) IP1 = NI
	   IM1  = I-1
	   IF(i.EQ.1) IM1 =1
         do j=2,NJ-1
	    do k=2,NK-1

	   	AVG = aux*(PROP(IM1,j,k)+PROP(IP1,j,k)+PROP(i,j-1,k)+ &
                     PROP(i,j+1,k)+PROP(i,j,k-1)+PROP(i,j,k+1))

	   	STORE(i,j,k)=SFM1*PROP(i,j,k)+SF*AVG

             end do
          end do
       end do

!
!********************************************************************************
!  Now smooth the property on the upper and lower boundaries.
!*******************************************************************************
! Considering the points in Y-frontier:

         do i=2,NI-1
	    do k=2,NK-1

      		AVG1  = aux2*(PROP(i-1,1,k)+PROP(i+1,1,k)+PROP(i,1,k+1)+PROP(i,1,k-1)+ &
                         2.*PROP(i,2,k) - PROP(i,3,k))

      		AVGNJ = aux2*(PROP(i-1,NJ,k)+PROP(i+1,NJ,k)+PROP(i,NJ,k+1)+PROP(i,NJ,k-1)+ &
                         2.*PROP(i,NJ-1,k) - PROP(i,NJ-2,k))

		STORE(i,1,k)=SFM1*PROP(i,1,k)+SF*AVG1

		STORE(i,NJ,k)=SFM1*PROP(i,NJ,k)+SF*AVGNJ

	    end do
         end do

! Considering the points in Z-frontier:

         do i=2,NI-1
	    do j=2,NJ-1

      		AVG2  = aux2*(PROP(i-1,j,1)+PROP(i+1,j,1)+PROP(i,j+1,1)+PROP(i,j-1,1)+ &
                        2.*PROP(i,j,2) - PROP(i,j,3))

      		AVGNK = aux2*(PROP(i-1,j,NK)+PROP(i+1,j,NK)+PROP(i,j+1,NK)+PROP(i,j-1,NK)+ &
                         2.*PROP(i,j,NK-1) - PROP(i,j,NK-2))

		STORE(i,j,1)=SFM1*PROP(i,j,1)+SF*AVG2

		STORE(i,j,NK)=SFM1*PROP(i,j,NK)+SF*AVGNK

	    end do
         end do
!
!----------------------------------------------------------------------------------
!
! Considering the edges in X-direction:

         do i=1,NI
	   IP1  = i+1
	   IF(i.EQ.NI) IP1 = NI
	   IM1  = I-1
	   IF(i.EQ.1) IM1 =1

	   	AVG = aux3*(PROP(IP1,1,1)+PROP(IM1,1,1)+ 2*PROP(i,2,1)- &
                     PROP(i,3,1)+ 2*PROP(i,1,2)-PROP(i,1,3))
	   	STORE(i,1,1)=SFM1*PROP(i,1,1)+SF*AVG
! -------------
	   	AVG = aux3*(PROP(IP1,NJ,1)+PROP(IM1,NJ,1)+ 2*PROP(i,NJ-1,1)- &
                     PROP(i,NJ-2,1)+ 2*PROP(i,NJ,2)-PROP(i,NJ,3))
	   	STORE(i,NJ,1)=SFM1*PROP(i,NJ,1)+SF*AVG
!--------------
	   	AVG = aux3*(PROP(IP1,1,NK)+PROP(IM1,1,NK)+ 2*PROP(i,2,NK)- &
                     PROP(i,3,NK)+ 2*PROP(i,1,NK-1)-PROP(i,1,NK-2))
	   	STORE(i,1,NK)=SFM1*PROP(i,1,NK)+SF*AVG
!--------------
	   	AVG = aux3*(PROP(IP1,NJ,NK)+PROP(IM1,NJ,NK)+ 2*PROP(i,NJ-1,NK)- &
                     PROP(i,NJ-2,NK)+ 2*PROP(i,NJ,NK-1)-PROP(i,NJ,NK-2))
	   	STORE(i,NJ,NK)=SFM1*PROP(i,NJ,NK)+SF*AVG

         end do
!
! Considering the edges in Y-direction:

         do j=2,NJ-1
	   	AVG1 = aux3*(PROP(1,j+1,1)+PROP(1,j-1,1)+ 2*PROP(2,j,1)- &
                     PROP(3,j,1)+ 2*PROP(1,j,2)-PROP(1,j,3))
	   	STORE(1,j,1)=SFM1*PROP(1,j,1)+SF*AVG1
! -------------
	   	AVG1 = aux3*(PROP(NI,j+1,1)+PROP(NI,j-1,1)+ 2*PROP(NI-1,j,1)- &
                     PROP(NI-2,j,1)+ 2*PROP(NI,j,2)-PROP(NI,j,3))
	   	STORE(NI,j,1)=SFM1*PROP(NI,j,1)+SF*AVG1
!--------------
	   	AVG1 = aux3*(PROP(1,j+1,NK)+PROP(1,j-1,NK)+ 2*PROP(2,j,NK)- &
                     PROP(3,j,NK)+ 2*PROP(1,j,NK-1)-PROP(1,j,NK-2))
	   	STORE(1,j,NK)=SFM1*PROP(1,j,NK)+SF*AVG1
!--------------
	   	AVG1 = aux3*(PROP(NI,j+1,NK)+PROP(NI,j-1,NK)+ 2*PROP(NI-1,j,NK)- &
                     PROP(NI-2,j,NK)+ 2*PROP(NI,j,NK-1)-PROP(NI,j,NK-2))
	   	STORE(NI,j,NK)=SFM1*PROP(NI,j,NK)+SF*AVG1

         end do
!
! Considering the edges in Z-direction:

         do k=2,NK-1
	   	AVG2 = aux3*(PROP(1,1,k+1)+PROP(1,1,k-1)+ 2*PROP(2,1,k)- &
                     PROP(3,1,k)+ 2*PROP(1,2,k)-PROP(1,3,k))
	   	STORE(1,1,k)=SFM1*PROP(1,1,k)+SF*AVG2
! -------------
	   	AVG2 = aux3*(PROP(NI,1,k+1)+PROP(NI,1,k-1)+ 2*PROP(NI-1,1,k)- &
                     PROP(NI-2,1,k)+ 2*PROP(NI,2,k)-PROP(NI,3,k))
	   	STORE(NI,1,k)=SFM1*PROP(NI,1,k)+SF*AVG2
!--------------
	   	AVG2 = aux3*(PROP(1,NJ,k+1)+PROP(1,NJ,k-1)+ 2*PROP(2,NJ,k)- &
                     PROP(3,NJ,k)+ 2*PROP(1,NJ-1,k)-PROP(1,NJ-2,k))
	   	STORE(1,NJ,k)=SFM1*PROP(1,NJ,k)+SF*AVG2
!--------------
	   	AVG2 = aux3*(PROP(NI,NJ,k+1)+PROP(NI,NJ,k-1)+ 2*PROP(NI-1,NJ,k)- &
                     PROP(NI-2,NJ,k)+ 2*PROP(NI,NJ-1,k)-PROP(NI,NJ-2,k))
	   	STORE(NI,NJ,k)=SFM1*PROP(NI,NJ,k)+SF*AVG2
         end do
!
! Reset the smoothed value to "PROP" before returning to the main program.

       do i=2,NI-1
          do j=2,NJ-1
	     do k=2,NK-1
                  IF (POR_VOL(i,j,k).NE.0.0) THEN

	   	     PROP(i,j,k) = STORE(i,j,k)

                 END IF
             end do
          end do
       end do

       RETURN
       END
