subroutine set_timestep

       use common_block

       implicit none
!
       real (kind = dp ) :: UMAX
!
!  This subroutine sets the length of the time step based on the
!  stagnation speed of sound "ASTAG" and the minimum length scale
!  of any element, "DMIN". The timestep must be called "DELTAT"

!  An assumption that the maximum flow speed will be equal to "ASTAG"
!  is also made. This will be pessimistic for subsonic flows
!  but may be optimistic for supersonic flows. In the latter case the
!  length of the time step as determined by "CFL" may need to be reduced.
!
!  The CFL number was input as data in data set "FLOW"

       ASTAG  = SQRT(GAMMA*RGAS*TSTAGIN)
       UMAX   = ASTAG

!   Calculate "deltat".

       DELTAT = CFL*DMIN/(UMAX+ASTAG)


       open(unit=35,file='TIME_STEP.dat')
       write(35,"(A)")('UMAX CFL DMIN DELTAT')
       write(35,"(4F10.5)") UMAX, CFL, DMIN, DELTAT

       RETURN
       END
