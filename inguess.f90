subroutine inguess

!  You do not need to touch this subroutine.

       use common_block

       implicit none

       INTEGER i,j,k


!  Set some completely arbitary values to enough of the flow
!  variables to enable plotting of the grid. These are not used at all
!  in the later calculations.

       do i=1,NI
          do j=1,NJ
             do k=1,NK
	   	RO(i,j,k)    = 1.2
	   	ROVX(i,j,k)  = 100.*FLOAT(i)/NI
	   	ROVY(i,j,k)  = 0.0
	   	ROVZ(i,j,k)  = 0.0
	   	P(i,j,k)     = 100000.*(0.9 + 0.1*FLOAT(i)/NI)
	   	HSTAG(i,j,k) = 300000.
	   	ROE(i,j,k)   = P(i,j,k)/(GAMMA-1.)
             end do
          end do
       end do

!     CALCULATE THE REFERENCE VALUES WHICH ARE USED TO CHECK CONVERGENCE

       ROIN      = PSTAGIN/RGAS/TSTAGIN
       REF_RO    = (PSTAGIN-PDOWN)/RGAS/TSTAGIN
       REF_T     = TSTAGIN*(PDOWN/PSTAGIN)**FGA
       REF_V     = SQRT(2*CP*(TSTAGIN-REF_T))
       REF_ROVX   = ROIN*REF_V
       REF_ROVY   = REF_ROVX
       REF_ROVZ   = REF_ROVX
       REF_ROE   = ROIN*CV*(TSTAGIN-REF_T)

       RETURN
       END
