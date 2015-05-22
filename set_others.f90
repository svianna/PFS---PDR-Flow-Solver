subroutine set_others

!  This routine calculates secondary flow variables from the primary ones
!  at every grid point.

       use common_block

       implicit none
!
       integer ( kind = 4 ) :: i,j,k
       real (kind = dp ) :: dummy,TEMP
!
!   The primary variables are RO,ROVX,ROVY and ROE

!  The secondary variables are the velocity components VX(i,j,k) and VY(i,j,k),
!  the static pressure P(i,j,k) and the stagnation enthalpy HSTAG(i,j,k).
!  Note:  "HSTAG"  NOT  "HO".

       do i=1,NI
          do j=1,NJ
             do k=1,NK     
			
	     	VX(i,j,k)=ROVX(i,j,k)/RO(i,j,k)
	    	VY(i,j,k)=ROVY(i,j,k)/RO(i,j,k)
	     	VZ(i,j,k)=ROVZ(i,j,k)/RO(i,j,k)
             	dummy=(VX(i,j,k)**2.+VY(i,j,k)**2.+VZ(i,j,k)**2.)
	     	TEMP=((ROE(i,j,k)/RO(i,j,k))-0.5*dummy)/CV
	     	P(i,j,k)=RO(i,j,k)*RGAS*TEMP
	     	HSTAG(i,j,k)=(CP*TEMP)+(0.5*dummy)

             end do
          end do
       end do
       

       RETURN
       END
