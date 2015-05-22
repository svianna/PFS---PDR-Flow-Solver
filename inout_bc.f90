subroutine inout_bc

!  This subroutine applies the boundary conditions that P = PDOWN
!  at I = NI. At the inlet boundary the change in density is relaxed
!  to obtain ROINLET(J) which is then used to obtain the other properties
!  at inlet assuming isentropic flow from stagnation conditions "PSTAGIN"
!  and TSTAGIN" together with the specified inlet flow angle "ALPHA1".

       use common_block

       implicit none

       integer ( kind = 4 ) :: i,j,k
       real (kind = dp ) :: aux, aux1, RFIN, RFIN1, ROSTAGIN

!  Because the inlet condition may become unstable it is safer to
!  relax the changes in inlet density by a factor "RFIN"
!  Typically "RFIN" = 0.25 as set below. Reduce this if the inlet
!  becomes unstable.
!
!  It is also worth checking if "ROINLET" is greater than "ROSTAGIN"
!  and setting ROINLET to 0.9999*ROSTAGIN if it is.
!  This saves the program crashing during transients.


       RFIN     = 0.25
       RFIN1    = 1.0-RFIN
       ROSTAGIN = PSTAGIN/(RGAS*TSTAGIN)


!********************************************************************************
! ---> Considering the faces in X-direction:
!
         do j=1,NJ
           do k=1,NK

	      IF(NSTEP.EQ.1) THEN
	  	 ROINLET(j,k) = RO(1,j,k)
	      ELSE
	         ROINLET(j,k) = RFIN*RO(1,j,k) + RFIN1*ROINLET(j,k)
              ENDIF

              IF (ROINLET(j,k).GT.0.9999*ROSTAGIN) ROINLET(j,k)=0.9999*ROSTAGIN
      
!      INSERT YOUR CODE HERE to calculate P(1,J),ROVX(1,J),ROVY(1,J)
!      and ROE(1,J)  from ROINLET(J), PSTAGIN, TSTAGIN  and ALPHA1.
!      also set VX(1,J), VY(1,J) and HSTAG(1,J)
      
	      aux=(ROSTAGIN/ROINLET(j,k))**GAMMA
   
	      P(1,j,k)=PSTAGIN/aux
!	      aux1=(PSTAGIN/P(1,J))**FGA
	      T(1,j,k)=P(1,j,k)/ROINLET(j,k)/RGAS
!	      V(1,j,k)=sqrt(2.*CP*(TSTAGIN-T(1,j,k)))

	      V(1,j,k)=sqrt(2*(PSTAGIN-P(1,j,k)/ROINLET(j,k)))
	      V(1,j,k)=sqrt(2.*cp*(TSTAGIN-T(1,j,k)))
	      VX(1,j,k)=V(1,j,k)
	      VY(1,j,k)=0.0
	      VZ(1,j,k)=0.0
	      ROE(1,j,k)=ROINLET(j,k)*(CV*T(1,j,k)+0.5*V(1,j,k)**2.)	
	      ROVX(1,j,k)=ROINLET(j,k)*VX(1,j,k)
	      ROVY(1,j,k)=ROINLET(j,k)*VY(1,j,k)
	      ROVZ(1,j,k)=ROINLET(j,k)*VZ(1,j,k)
	      HSTAG(1,j,k)=(CP*T(1,j,k))+(0.5*V(1,j,k)*V(1,j,k))

           end do
         end do
!
!---------------------------------------------------------------------------------
!
!**********************************************************************************
!**********************************************************************************
!   Set the pressure at the downstream boundaries to the exit
!   static pressure "PDOWN" for all j values.
!
          do j=1,NJ
	     do k=1,NK
		P(NI,j,k)=PDOWN
	     end do
          end do
!
! ------------------------------------------------------------

       RETURN
       END
