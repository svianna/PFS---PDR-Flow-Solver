subroutine output(DUMP)

! This subroutine writes a file to unit 7 for use by the plotting program "EULPLT".

       use common_block

       implicit none

       CHARACTER FILENAME*16, FILENAME2*16, NAMEDUMP
       integer ( kind = 4 ) :: i,j,k,IMAX,JMAX,KMAX,DUMP,ID,IDVTK
       integer ( kind = 4 ) :: NCELLSS,count,countt
       integer ( kind = 4 ) :: COUNT1,COUNT2,COUNT3,COUNT4,COUNT5
       integer ( kind = 4 ) :: aux2, aux3,aux4,aux5
       real (kind = dp ) :: EKE,TSTAT,VELSQ,TSTAG
       real (kind = dp ) :: PSTAG(NI,NJ,NK), VMACH(NI,NJ,NK)
!
!
!     CALCULATE THE SECONDARY VARIABLES
!
       do i=1,NI
          do j=1,NJ
             do k=1,NK

		VX(i,j,k)=ROVX(i,j,k)/RO(i,j,k)
		VY(i,j,k)=ROVY(i,j,k)/RO(i,j,k)
		VZ(i,j,k)=ROVZ(i,j,k)/RO(i,j,k)
		EKE=0.5*(VX(i,j,k)**2.+VY(i,j,k)**2.+VZ(i,j,k)**2.)
		TSTAT=(HSTAG(i,j,k) - EKE)/CP
		P(i,j,k)=RO(i,j,k)*RGAS*TSTAT
		ROE(i,j,k) = RO(i,j,k)*(CV*TSTAT + EKE)

             end do
          end do
       end do
!
!    CALCULATE THE MACH NUMBER AND STAGNATION PRESSSURE
!
       do i=1,NI
          do j=1,NJ
             do k=1,NK

		TSTAT = P(i,j,k)/RGAS/RO(i,j,k)
		VELSQ = VX(i,j,k)**2.+VY(i,j,k)**2.+VZ(i,j,k)**2.
		TSTAG = TSTAT + 0.5*VELSQ/CP
		VMACH(i,j,k) = SQRT(VELSQ/(GAMMA*RGAS*TSTAT))
		PSTAG(i,j,k) = P(i,j,k)*(TSTAG/TSTAT)**(1/FGA)

             end do
          end do
       end do
!
!************************************************************************* 

        IDVTK=21

          WRITE(FILENAME,600)NSTEP

 600    FORMAT('D_',I7,'.vtk')
        DO I=1,11
          IF(FILENAME(I:I).EQ.' ')THEN
            FILENAME(I:I)='0'
          ENDIF
        ENDDO              
        
        OPEN(UNIT=IDVTK,FILE=FILENAME)


       write(IDVTK,'(A,3I4)')'# vtk DataFile Version 3.0'
       write(IDVTK,'(A)') 'Mesh Files'
       write(IDVTK,'(A,3I4)') 'ASCII'
       write(IDVTK,'(A,3I4)') 'DATASET UNSTRUCTURED_GRID'
       write(IDVTK,'(A,1I10,A)') 'POINTS', NPOINTS,   ' float'
          NCELLSS= (NI-1)*(NJ-1)*(NK-1) 
          aux3=NCELLSS*9.
          DO i=1,NI
             DO j=1,NJ
                DO k=1,NK
                     write(IDVTK,'(3F7.2)') X(i,j,k),Y(i,j,k),Z(i,j,k)
                END DO
             END DO
          END DO
          write(IDVTK,*)

          write(IDVTK,'(A,2I10)')'CELLS',NCELLSS,aux3

! escrevendo a connectividade

       aux2 = 0
       aux3 = aux2 + 1
       aux4 = NJ*NK
       aux5 =(NJ-1)*(NK-1) 
     
      DO i=1, NCELLSS

         WRITE(IDVTK,'(9I6)') 8,aux2,aux4+aux2,aux4+aux2+NK,aux2+NK, &
                           aux2+1,aux4+aux2+1,aux4+NK+aux2+1,aux2+NK+1

         aux2 = aux2 + 1.
         aux3 = aux3 + 1.

       	if (MOD(i,(NK-1)).eq.0) then
           aux3 = aux3 + 1.
           aux2 = aux2 + 1.
       	endif

       	if (MOD(i,aux5).eq.0) then
           aux3 = aux3 + NK
           aux2 = aux2 + NK
       	endif
      
      END DO
      
         write(IDVTK, *) 
         write(IDVTK,'(A,1I10)')'CELL_TYPES',NCELLSS
      	 do i=1,NCELLSS
             WRITE(21,'(1I2)') 12
      	 end do

       write(IDVTK,*) 'POINT_DATA', NPOINTS
       write (IDVTK,'(A,1I2)') 'SCALARS DENSITY float', 1
       write (IDVTK,'(A)') 'LOOKUP_TABLE default'
          
       do i=1,NI
          do j=1,NJ
             do k=1,NK
                write (IDVTK,201) RO(i,j,k) 
	     end do
         end do
       end do
!****************************************************************  
            
       write(IDVTK,'(A)') 'VECTORS velocity_x float'
       do i=1,NI
          do j=1,NJ
             do k=1,NK
               write (IDVTK,201) VX(i,j,k),0.0, 0.0
	     end do
          end do
       end do

       write(IDVTK,'(A)') 'VECTORS velocity_y float'
       do i=1,NI
          do j=1,NJ
             do k=1,NK
               write (IDVTK,201)  0.0, VY(i,j,k), 0.0
             end do
          end do
       end do

       write(IDVTK,'(A)') 'VECTORS velocity_z float'
       do i=1,NI
          do j=1,NJ
             do k=1,NK
               write (IDVTK,201) 0.0, 0.0, VZ(i,j,k)
             end do
          end do
       end do
         
      write (IDVTK,'(A,1I2)') 'SCALARS PRESSURE float', 1
      write (IDVTK,'(A)') 'LOOKUP_TABLE default'
       do i=1,NI
          do j=1,NJ
             do k=1,NK
                write(IDVTK,202) P(i,j,k) 
              end do
           end do
        end do
      
       write (IDVTK,'(A,1I2)') 'SCALARS MACH_NUMBER float', 1
       write (IDVTK,'(A)') 'LOOKUP_TABLE default'
       do i=1,NI
          do j=1,NJ
             do k=1,NK
                write(IDVTK,201) VMACH(i,j,k) 
              end do
           end do
        end do

!      WRITE(21,'(A)') 'VECTORS DEL_ROVX float'
!          do k=1,(NK+1)
!            do i=1,(NI+1)
!               do j=1,(NJ+1)
!               WRITE(21,203) DELROVX(i,j,k), 0.0, 0.0, 0.0
!               end do
!            end do
!         end do

       close(IDVTK)

!*************************************************************************


     ! write(20,*)  'X VELOCITY'    
      !write(20,201)(((VX(i,j,k),i=1,DNI),j=1,1),k=1,1)

  201 FORMAT(10F10.4)

      !write(20,*)  ' Y  VELOCITY'     
      !write(20,201)(((VY(i,j,k),i=1,1),j=1,DNJ),k=1,1)

      !write(20,*)  ' Z  VELOCITY'     
      !write(20,201)(((VZ(i,j,k),i=1,1),j=1,1),k=1,DNK)

      !write(20,*)  ' MACH NUMBER '    
      !write(20,201)(((VMACH(i,j,k),i=1,DNI),j=1,DNJ),k=1,DNK)

      !write(20,*)  ' STATIC PRESSURE'   
      !write(20,202)(((P(i,j,k),i=1,DNI),j=1,DNJ),k=1,DNK)

  202 FORMAT(10F10.1)

      !write(20,*)  ' DENSITY '   
      !write(20,201)(((RO(i,j,k),i=1,DNI),j=1,DNJ),k=1,DNK)

!      write(20,*)  ' VARIABLE PLOTVAR '   
!      write(20,201)(((PLOTVAR(i,j,k),i=1,DNI),j=1,DNJ),k=1,DNK)

!      write(20,*)  ' DEL ROVX '   
!      write(20,203)(((DELROVX(i,j,k),i=1,DNI),j=1,DNJ),k=1,DNK)
!  203 FORMAT(10E10.4)

       CLOSE(20)

       RETURN
       END
