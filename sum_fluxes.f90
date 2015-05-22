subroutine sum_fluxes(IFLUX,JFLUX,KFLUX,PROP,DELPROP)

!  This subroutine sums the fluxes for each element, calculates
!  the changes in the variable "PROP" and distributes them to the
!  four corners of the element.

       use common_block
!
       implicit none
!
       integer ( kind = 4 ) :: i,j,k
       real (kind = dp ) :: ADD
       real (kind = dp ) :: IFLUX(NI,(NJ-1),(NK-1)),JFLUX((NI-1),NJ,(NK-1)), &
                             KFLUX((NI-1),(NJ-1),NK), STORE((NI-1),(NJ-1),(NK-1)),&
                             PROP(NI,NJ,NK),DELPROP(NI,NJ,NK)
                

!
!   Find the change in the variable "PROP" in each cell over the
!   time step "DELTAT" and save it in "STORE(I,J)"
!

      DO I=1,NI-1
      DO J=1,NJ-1
      DO K=1,NK-1

        IF (POR_VOL(i,j,k).NE.0.0) THEN
!       
           store(i,j,k)=(iflux(i,j,k)-iflux(i+1,j,k)+ jflux(i,j,k)- &
           jflux(i,j+1,k)+ kflux(i,j,k)-kflux(i,j,k+1))*(deltat/  &
             volcell(i,j,k)*POR_VOL(i,j,k))

        ELSE
!       
          store(i,j,k)=(iflux(i,j,k)-iflux(i+1,j,k)+ jflux(i,j,k)- &
          jflux(i,j+1,k)+ kflux(i,j,k)-kflux(i,j,k+1))*(deltat/volcell(i,j,k))

        END IF   
!
!      calculate the change in the variable
!     "PROP" over the time step "deltat" and sets the result to "STORE(I,J)"
!
      END DO
      END DO
      END DO
!

!  Now distribute the changes equally to the four interior corners of each
!  cell. Each interior grid points receive one quarter of the change
!  from each of the four cells adjacent to it.
!
      DO I=2,NI-1
      DO J=2,NJ-1
      DO K=2,NK-1

      
         ADD=0.125*(store(i,j-1,k)+store(i-1,j-1,k)+store(i,j-1,k-1)+store(i-1,j-1,k-1)+ &
                    store(i,j,k)+store(i-1,j,k)+store(i,j,k-1)+store(i-1,j,k-1))
 !
        IF (POR_VOL(i,j,k).NE.0.0) THEN
!       
         PROP(I,j,K) = PROP(I,j,K)     + ADD

        END IF

     END DO
     END DO
     END DO
!
!--------------------------------------------------------------------------------------
!  Now add the changes to the boundaries.
!  These receive one quarter of the change from each of the four cells
!  adjacent to them.
!
!On Y faces
!
      DO I=2,NI-1
      DO K=2,NK-1
      
      ADD=0.25*(store(i,1,k)+store(i-1,1,k)+store(i,1,k-1)+store(i-1,1,k-1))         
      PROP(I,1,K)     = PROP(I,1,K)     + ADD
!
      ADD=0.25*(store(i,NJ-1,k)+store(i-1,NJ-1,k)+store(i,NJ-1,k-1)+store(i-1,NJ-1,k-1))         
      PROP(I,NJ,K)     = PROP(I,NJ,K)     + ADD

      END DO
      END DO
!
!On Z faces
!
      DO I=2,NI-1
      DO J=2,NJ-1
      
      ADD=0.25*(store(i,j,1)+store(i-1,j,1)+store(i,j-1,1)+store(i-1,j-1,1))         
      PROP(I,J,1)     = PROP(I,J,1)     + ADD
!
      ADD=0.25*(store(i,j,NK-1)+store(i-1,j,NK-1)+store(i,j-1,NK-1)+store(i-1,j-1,NK-1))         
      PROP(I,J,NK)     = PROP(I,J,NK)     + ADD

      END DO
      END DO
!
!On X faces
!
      DO J=2,NJ-1
      DO K=2,NK-1
      
      ADD=0.25*(store(1,j,k)+store(1,j-1,k)+store(1,j,k-1)+store(1,j-1,k-1))         
      PROP(1,J,K)     = PROP(1,J,K)     + ADD
!
      ADD=0.25*(store(NI-1,j,k)+store(NI-1,j-1,k)+store(NI-1,j,k-1)+store(NI-1,j-1,k-1))         
      PROP(NI,J,K)     = PROP(NI,J,K)     + ADD

      END DO
      END DO
!
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!  Now add the changes to the boundary edges.
!  These receive half the change from each of the two cells
!  adjacent to them.
!
!On X direction
!
      DO I=2,NI-1
!
      ADD=0.5*(STORE(I,1,1)+STORE(I-1,1,1))
      PROP(I,1,1)     = PROP(I,1,1)    +  ADD
!
      ADD=0.5*(STORE(I,1,NK-1)+STORE(I-1,1,NK-1))
      PROP(I,1,NK)    = PROP(I,1,NK)   +  ADD
!
      ADD=0.5*(STORE(I,NJ-1,1)+STORE(I-1,NJ-1,1))
      PROP(I,NJ,1)     = PROP(I,NJ,1)    +  ADD
!
      ADD=0.5*(STORE(I,NJ-1,NK-1)+STORE(I-1,NJ-1,NK-1))
      PROP(I,NJ,NK)    = PROP(I,NJ,NK)   +  ADD
!
      END DO
!
!On Z direction
!
      DO K=2,NK-1
!
      ADD=0.5*(STORE(NI-1,1,K)+STORE(NI-1,1,K-1))
      PROP(NI,1,K)     = PROP(NI,1,K)   + ADD
!
      ADD=0.5*(STORE(1,1,K)+STORE(1,1,K-1))
      PROP(1,1,K)      = PROP(1,1,K)    + ADD
!
      ADD=0.5*(STORE(NI-1,NJ-1,K)+STORE(NI-1,NJ-1,K-1))
      PROP(NI,NJ,K)     = PROP(NI,NJ,K)   + ADD
!
      ADD=0.5*(STORE(1,NJ-1,K)+STORE(1,NJ-1,K-1))
      PROP(1,NJ,K)      = PROP(1,NJ,K)    + ADD
!
      END DO
!
!On Y direction
!
      DO J=2,NJ-1
!
      ADD=0.5*(STORE(NI-1,J,1)+STORE(NI-1,J-1,1))
      PROP(NI,J,1)     = PROP(NI,J,1)   + ADD
!
      ADD=0.5*(STORE(1,J,1)+STORE(1,J-1,1))
      PROP(1,J,1)      = PROP(1,J,1)    + ADD
!
      ADD=0.5*(STORE(NI-1,J,NK-1)+STORE(NI-1,J-1,NK-1))
      PROP(NI,J,NK)     = PROP(NI,J,NK)   + ADD
!
      ADD=0.5*(STORE(1,J,NK-1)+STORE(1,J-1,NK-1))
      PROP(1,J,NK)      = PROP(1,J,NK)    + ADD
!
      END DO
!
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!      Now add the changes on to the eight corner points.
!      These receive the full change from the single cell of which
!      they form one corner.
!
      ADD=STORE(1,1,1)
      PROP(1,1,1)       = PROP(1,1,1)  + ADD
!
      ADD=STORE(1,1,NK-1)
      PROP(1,1,NK)      = PROP(1,1,NK)  + ADD 
!      
      ADD=STORE(NI-1,1,1)
      PROP(NI,1,1)      = PROP(NI,1,1) + ADD
!
      ADD=STORE(NI-1,1,NK-1)
      PROP(NI,1,NK)     = PROP(NI,1,NK) + ADD
!
      ADD=STORE(1,NJ-1,1)
      PROP(1,NJ,1)       = PROP(1,NJ,1)  + ADD
!
      ADD=STORE(1,NJ-1,NK-1)
      PROP(1,NJ,NK)      = PROP(1,NJ,NK)  + ADD 
!      
      ADD=STORE(NI-1,NJ-1,1)
      PROP(NI,NJ,1)      = PROP(NI,NJ,1) + ADD
!
      ADD=STORE(NI-1,NJ-1,NK-1)
      PROP(NI,NJ,NK)     = PROP(NI,NJ,NK) + ADD
!
! Now save the changes in the primary variables as "DELPROP".
! This will be used in the convergence check and also in future
! improvements to the scheme.
!
      DO  I=1,NI-1
      DO  J=1,NJ-1
      DO  K=1,NK-1
      DELPROP(I,J,K) = STORE(I,J,K)
      END DO
      END DO
      END DO
!
!
      RETURN
      END
