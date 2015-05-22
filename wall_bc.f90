subroutine wall_bc

! This subroutine can read the data on the mesh and flow conditions.

       use common_block

       implicit none
!
       integer ( kind = 4 ) :: i,j,k
!
!***********************************************************************
!       Correcting the fluxes at the boundaries
!***********************************************************************
!
!       print*,'WALL_BC'
! ----------------------------------------------------
!
          do i=1,NI-1
	     do k=1,NK-1
		FLUXJ_MA(i,1,k)=0.0
		FLUXJ_MA(i,NJ,k)=0.0
	     end do
          end do
!
! ----------------------------------------------------
!
          do i=1,NI-1
	     do j=1,NJ-1
		FLUXK_MA(i,j,1)=0.0
		FLUXK_MA(i,j,NK)=0.0
	     end do
          end do
!
!*************************************************************************
!!-----------------------------------------------------------------------
!!
!!TESTE
       do i=1,NI-1
	  do j=1,NJ-1
	      do k=1,NK-1
            
                 if (POR_VOL(i,j,k).ne.1) then
                     FLUXI_MA(i,j,k)=0.0
                     FLUXJ_MA(i,j,k)=0.0
                     FLUXK_MA(i,j,k)=0.0
                 end if

	       end do
	  end do
       end do

!------------------------------------------------------------------------
!
!
       return
       end
