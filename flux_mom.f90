subroutine  flux_mom

!  This subroutine calculates the fluxes of momentum
!  across the faces of every cell.

!  The faces with i,j and k = 1 or ni,nj and nk are the upstream and downstream
!  boundaries to the flow domain. All fluxes are calculated assuming a linear
!  variation in the flow properties between the cell corner nodes.

       use common_block

       implicit none
!
       integer ( kind = 4 ) :: i,j,k
       real (kind = dp ) :: res1, res2, res3, res4, res
       real (kind = dp ) :: AUX
!
! *******************************************
!       RECALCULATE THE FLUX OF MASS
! *******************************************
!
       do i=1,NI
          do j=1,NJ-1
	     do k=1,NK-1

		FLUXI_MA(i,j,k)=0.25*( (ROVX(i,j,k)+ROVX(i,j,k+1)+ROVX(i,j+1,k)+ &
                         ROVX(i,j+1,k+1))*DAIX(i,j,k)*POR_AX(i,j,k) + (ROVY(i,j,k)+ROVY(i,j,k+1)+ &
                         ROVY(i,j+1,k)+ROVY(i,j+1,k+1))*DAIY(i,j,k)*POR_AX(i,j,k) + (ROVZ(i,j,k)+ &
                         ROVZ(i,j,k+1)+ROVZ(i,j+1,k)+ROVZ(i,j+1,k+1))*DAIZ(i,j,k)*POR_AX(i,j,k) )

             end do
          end do
       end do
!
!
       do i=1,NI-1
          do j=1,NJ
	     do k=1,NK-1

		FLUXJ_MA(i,j,k)=0.25*( (ROVX(i,j,k)+ROVX(i+1,j,k)+ROVX(i+1,j,k+1)+ &
                         ROVX(i,j,k+1))*DAJX(i,j,k)*POR_AY(i,j,k) + (ROVY(i,j,k)+ROVY(i+1,j,k)+ &
                         ROVY(i+1,j,k+1)+ROVY(i,j,k+1))*DAJY(i,j,k)*POR_AY(i,j,k) + (ROVZ(i,j,k)+ &
                         ROVZ(i+1,j,k)+ROVZ(i+1,j,k+1)+ROVZ(i,j,k+1))*DAJZ(i,j,k)*POR_AY(i,j,k) )
             end do
          end do
       end do


       do i=1,NI-1
          do j=1,NJ-1
	     do k=1,NK

		FLUXK_MA(i,j,k)=0.25*( (ROVX(i,j,k)+ROVX(i+1,j,k)+ROVX(i+1,j+1,k)+ &
                         ROVX(i,j+1,k))*DAKX(i,j,k)*POR_AZ(i,j,k) + (ROVY(i,j,k)+ROVY(i+1,j,k)+ &
                         ROVY(i+1,j+1,k)+ROVY(i,j+1,k))*DAKY(i,j,k)*POR_AZ(i,j,k) + (ROVZ(i,j,k)+ &
                         ROVZ(i+1,j,k)+ROVZ(i+1,j+1,k)+ROVZ(i,j+1,k))*DAKZ(i,j,k)*POR_AZ(i,j,k) )     
             end do
          end do
       end do

! ----------------------------------------------------
!       Correcting the fluxes at the boundaries
! ----------------------------------------------------
!
      call wall_bc
!
! *******************************************
!       CALCULATE THE FLUX OF MOMENTUM
! *******************************************

!	CD = 1.3
        AUX = - 0.5*CD


! ---> IN X DIRECTION --> X-MOMENTUM

! Component i
       do i=1,NI 
          do j=1,NJ-1
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AX(i,j,k))*(DAIX(i,j,k))*RO(i,j,k)*ABS(VX(i,j,k))*VX(i,j,k)
        	res2 = AUX*(1-POR_AX(i,j,k))*(DAIX(i,j,k))*RO(i,j+1,k)*ABS(VX(i,j+1,k))*VX(i,j+1,k)
        	res3 = AUX*(1-POR_AX(i,j,k))*(DAIX(i,j,k))*RO(i,j,k+1)*ABS(VX(i,j,k+1))*VX(i,j,k+1)
        	res4 = AUX*(1-POR_AX(i,j,k))*(DAIX(i,j,k))*RO(i,j+1,k+1)*ABS(VX(i,j+1,k+1))*VX(i,j+1,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)

        	FLUXI_XMOM(i,j,k)=0.25*( FLUXI_MA(i,j,k)*(VX(i,j,k)+ &
                        VX(i,j,k+1)+VX(i,j+1,k)+VX(i,j+1,k+1))+ &
                        (P(i,j,k)+P(i,j,k+1)+P(i,j+1,k)+P(i,j+1,k+1))* &
                        DAIX(i,j,k)*POR_AX(i,j,k) ) + res

             end do
          end do
       end do

! Component j
       do i=1,NI-1 
          do j=1,NJ
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AY(i,j,k))*(DAJX(i,j,k))*RO(i,j,k)*ABS(VX(i,j,k))*VX(i,j,k)
        	res2 = AUX*(1-POR_AY(i,j,k))*(DAJX(i,j,k))*RO(i+1,j,k)*ABS(VX(i+1,j,k))*VX(i+1,j,k)
        	res3 = AUX*(1-POR_AY(i,j,k))*(DAJX(i,j,k))*RO(i,j,k+1)*ABS(VX(i,j,k+1))*VX(i,j,k+1)
        	res4 = AUX*(1-POR_AY(i,j,k))*(DAJX(i,j,k))*RO(i+1,j,k+1)*ABS(VX(i+1,j,k+1))*VX(i+1,j,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)


        	FLUXJ_XMOM(i,j,k) = 0.25*(FLUXJ_MA(i,j,k)*(VX(i,j,k)+ &
                        VX(i,j,k+1)+VX(i+1,j,k)+VX(i+1,j,k+1))+ &
                        (P(i,j,k)+P(i,j,k+1)+P(i+1,j,k)+P(i+1,j,k+1))* &
                        DAJX(i,j,k)*POR_AY(i,j,k) ) + res
             end do
          end do
       end do

! Component k
       do i=1,NI-1
          do j=1,NJ-1
	     do k=1,NK

        	res1 = AUX*(1-POR_AZ(i,j,k))*(DAKX(i,j,k))*RO(i,j,k)*ABS(VX(i,j,k))*VX(i,j,k)
        	res2 = AUX*(1-POR_AZ(i,j,k))*(DAKX(i,j,k))*RO(i,j+1,k)*ABS(VX(i,j+1,k))*VX(i,j+1,k)
        	res3 = AUX*(1-POR_AZ(i,j,k))*(DAKX(i,j,k))*RO(i+1,j,k)*ABS(VX(i+1,j,k))*VX(i+1,j,k)
        	res4 = AUX*(1-POR_AZ(i,j,k))*(DAKX(i,j,k))*RO(i+1,j+1,k)*ABS(VX(i+1,j+1,k))*VX(i+1,j+1,k)
        	res= 0.25*(res1 + res2 + res3 + res4)


       		FLUXK_XMOM(i,j,k)=0.25*( FLUXK_MA(i,j,k)*(VX(i,j,k)+ &
                        VX(i+1,j,k)+VX(i,j+1,k)+VX(i+1,j+1,k))+ &
                        (P(i,j,k)+P(i+1,j,k)+P(i,j+1,k)+P(i+1,j+1,k))* &
                        DAKX(i,j,k)*POR_AZ(i,j,k) ) + res

             end do
          end do
       end do


! ---> IN Y DIRECTION --> Y-MOMENTUM

! Component i
       do i=1,NI
          do j=1,NJ-1
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AX(i,j,k))*(DAIY(i,j,k))*RO(i,j,k)*ABS(VY(i,j,k))*VY(i,j,k)
        	res2 = AUX*(1-POR_AX(i,j,k))*(DAIY(i,j,k))*RO(i,j+1,k)*ABS(VY(i,j+1,k))*VY(i,j+1,k)
        	res3 = AUX*(1-POR_AX(i,j,k))*(DAIY(i,j,k))*RO(i,j,k+1)*ABS(VY(i,j,k+1))*VY(i,j,k+1)
        	res4 = AUX*(1-POR_AX(i,j,k))*(DAIY(i,j,k))*RO(i,j+1,k+1)*ABS(VY(i,j+1,k+1))*VY(i,j+1,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXI_YMOM(i,j,k)= 0.25*( FLUXI_MA(i,j,k)*(VY(i,j,k)+ &
                        VY(i,j+1,k)+VY(i,j+1,k+1)+VY(i,j,k+1))+ &
                        (P(i,j,k)+P(i,j+1,k)+P(i,j+1,k+1)+P(i,j,k+1))* &
                        DAIY(i,j,k)*POR_AX(i,j,k) ) + res
             end do
          end do
       end do

! Component j
       do i=1,NI-1
          do j=1,NJ
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AY(i,j,k))*(DAJY(i,j,k))*RO(i,j,k)*ABS(VY(i,j,k))*VY(i,j,k)
        	res2 = AUX*(1-POR_AY(i,j,k))*(DAJY(i,j,k))*RO(i+1,j,k)*ABS(VY(i+1,j,k))*VY(i+1,j,k)
        	res3 = AUX*(1-POR_AY(i,j,k))*(DAJY(i,j,k))*RO(i,j,k+1)*ABS(VY(i,j,k+1))*VY(i,j,k+1)
        	res4 = AUX*(1-POR_AY(i,j,k))*(DAJY(i,j,k))*RO(i+1,j,k+1)*ABS(VY(i+1,j,k+1))*VY(i+1,j,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXJ_YMOM(i,j,k)= 0.25*( FLUXJ_MA(i,j,k)*(VY(i,j,k)+ &
                        VY(i+1,j,k)+VY(i+1,j,k+1)+VY(i,j,k+1))+ &
                        (P(i,j,k)+P(i+1,j,k)+P(i+1,j,k+1)+P(i,j,k+1))* &
                        DAJY(i,j,k)*POR_AY(i,j,k) ) + res
             end do
          end do
       end do


! Component k
       do i=1,NI-1
          do j=1,NJ-1
	     do k=1,NK

        	res1 = AUX*(1-POR_AZ(i,j,k))*(DAKY(i,j,k))*RO(i,j,k)*ABS(VY(i,j,k))*VY(i,j,k)
        	res2 = AUX*(1-POR_AZ(i,j,k))*(DAKY(i,j,k))*RO(i+1,j,k)*ABS(VY(i+1,j,k))*VY(i+1,j,k)
        	res3 = AUX*(1-POR_AZ(i,j,k))*(DAKY(i,j,k))*RO(i,j+1,k)*ABS(VY(i,j+1,k))*VY(i,j+1,k)
        	res4 = AUX*(1-POR_AZ(i,j,k))*(DAKY(i,j,k))*RO(i+1,j+1,k)*ABS(VY(i+1,j+1,k))*VY(i+1,j+1,k)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXK_YMOM(i,j,k)= 0.25*( FLUXK_MA(i,j,k)*(VY(i,j,k)+ &
                        VY(i+1,j,k)+VY(i+1,j+1,k)+VY(i,j+1,k))+ &
                        (P(i,j,k)+P(i+1,j,k)+P(i+1,j+1,k)+P(i,j+1,k))* &
                        DAKY(i,j,k)*POR_AZ(i,j,k) ) + res
             end do
          end do
       end do

! ---> IN Z DIRECTION --> Z-MOMENTUM

! Component i
       do i=1,NI 
          do j=1,NJ-1
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AX(i,j,k))*(DAIZ(i,j,k))*RO(i,j,k)*ABS(VZ(i,j,k))*VZ(i,j,k)
        	res2 = AUX*(1-POR_AX(i,j,k))*(DAIZ(i,j,k))*RO(i,j,k+1)*ABS(VZ(i,j,k+1))*VZ(i,j,k+1)
        	res3 = AUX*(1-POR_AX(i,j,k))*(DAIZ(i,j,k))*RO(i,j+1,k)*ABS(VZ(i,j+1,k))*VZ(i,j+1,k)
        	res4 = AUX*(1-POR_AX(i,j,k))*(DAIZ(i,j,k))*RO(i,j+1,k+1)*ABS(VZ(i,j+1,k+1))*VZ(i,j+1,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXI_ZMOM(i,j,k)= 0.25*( FLUXI_MA(i,j,k)*(VZ(i,j,k)+ &
                        VZ(i,j,k+1)+VZ(i,j+1,k)+VZ(i,j+1,k+1))+ &
                        (P(i,j,k)+P(i,j,k+1)+P(i,j+1,k)+P(i,j+1,k+1))* &
                        DAIZ(i,j,k)*POR_AX(i,j,k) ) + res
	     end do
          end do
       end do

! Component j
       do i=1,NI-1
          do j=1,NJ
	     do k=1,NK-1

        	res1 = AUX*(1-POR_AY(i,j,k))*(DAJZ(i,j,k))*RO(i,j,k)*ABS(VZ(i,j,k))*VZ(i,j,k)
        	res2 = AUX*(1-POR_AY(i,j,k))*(DAJZ(i,j,k))*RO(i+1,j,k)*ABS(VZ(i+1,j,k))*VZ(i+1,j,k)
        	res3 = AUX*(1-POR_AY(i,j,k))*(DAJZ(i,j,k))*RO(i,j,k+1)*ABS(VZ(i,j,k+1))*VZ(i,j,k+1)
        	res4 = AUX*(1-POR_AY(i,j,k))*(DAJZ(i,j,k))*RO(i+1,j,k+1)*ABS(VZ(i+1,j,k+1))*VZ(i+1,j,k+1)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXJ_ZMOM(i,j,k)= 0.25*( FLUXJ_MA(i,j,k)*(VZ(i,j,k)+ &
                        VZ(i+1,j,k)+VZ(i,j,k+1)+VZ(i+1,j,k+1))+ &
                        (P(i,j,k)+P(i+1,j,k)+P(i,j,k+1)+P(i+1,j,k+1))* &
                        DAJZ(i,j,k)*POR_AY(i,j,k) ) + res
	     end do
          end do
       end do

! Component k
       do i=1,NI-1
          do j=1,NJ-1
	       do k=1,NK

        	res1 = AUX*(1-POR_AZ(i,j,k))*(DAKZ(i,j,k))*RO(i,j,k)*ABS(VZ(i,j,k))*VZ(i,j,k)
        	res2 = AUX*(1-POR_AZ(i,j,k))*(DAKZ(i,j,k))*RO(i+1,j,k)*ABS(VZ(i+1,j,k))*VZ(i+1,j,k)
        	res3 = AUX*(1-POR_AZ(i,j,k))*(DAKZ(i,j,k))*RO(i,j+1,k)*ABS(VZ(i,j+1,k))*VZ(i,j+1,k)
        	res4 = AUX*(1-POR_AZ(i,j,k))*(DAKZ(i,j,k))*RO(i+1,j+1,k)*ABS(VZ(i+1,j+1,k))*VZ(i+1,j+1,k)
        	res= 0.25*(res1 + res2 + res3 + res4)

		FLUXK_ZMOM(i,j,k)= 0.25*( FLUXK_MA(i,j,k)*(VZ(i,j,k)+ &
                        VZ(i+1,j,k)+VZ(i,j+1,k)+VZ(i+1,j+1,k))+ &
                        (P(i,j,k)+P(i+1,j,k)+P(i,j+1,k)+P(i+1,j+1,k))* &
                        DAKZ(i,j,k)*POR_AZ(i,j,k) ) + res

	     end do
          end do
       end do

       return
       END

