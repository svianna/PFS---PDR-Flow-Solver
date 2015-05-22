subroutine read_porosity

! This subroutine can read the porosity values.

!
       use common_block

       implicit none

	character (len = 64) :: output_file_name_1, output_file_name_2, output_file_name_3, &
							output_file_name_4, output_file_name_5, output_file_name_6, &
							output_file_name_7, str, word, tmp
	integer :: unit1, unit2, unit3, unit4, unit5, unit6, unit7

	integer, parameter :: sp = kind(0.0 )	! single precision
	integer :: i,j,k
!
!
       allocate ( POR_AX(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( POR_AY(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( POR_AZ(1:(NI-1),1:(NJ-1),1:NK) )
       allocate ( POR_VOL(1:(NI-1),1:(NJ-1),1:(NK-1)) )
!

	output_file_name_4 = 'cube.porosity.asc'
	unit4 = 4
	open (unit=unit4,file=output_file_name_4)
	do i = 1, DNI
		do j = 1, DNJ
			do k = 1, DNK
				read(unit4, '(f14.6)') POR_VOL(i,j,k)
			end do
		end do
	end do
	close (unit4)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_4


	output_file_name_1 = 'cube.porosity.areaX.asc'
	unit1 = 1
	open (unit=unit1,file=output_file_name_1)
	do i = 1, NI
		do j = 1, DNJ
			do k = 1, DNK
				read(unit1, '(f14.6)') POR_AX(i,j,k)
			end do
		end do
	end do
	close (unit1)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_1

	output_file_name_2 = 'cube.porosity.areaY.asc'
	unit2 = 2
	open (unit=unit2,file=output_file_name_2)
	do i = 1, DNI
		do j = 1, NJ
			do k = 1, DNK
				read(unit2, '(f14.6)') POR_AY(i,j,k)
			end do
		end do
	end do
	close (unit2)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_2

	output_file_name_3 = 'cube.porosity.areaZ.asc'
	unit3 = 3
	open (unit=unit3,file=output_file_name_3)
	do i = 1, DNI
		do j = 1, DNJ
			do k = 1, NK
				read(unit3, '(f14.6)') POR_AZ(i,j,k)
			end do
		end do
	end do
	close(unit3)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_3

	output_file_name_4 = 'cube.porosity.asc'
	unit4 = 4
	open (unit=unit4,file=output_file_name_4)
	do i = 1, DNI
		do j = 1, DNJ
			do k = 1, DNK
				read(unit4, '(f14.6)') POR_VOL(i,j,k)
			end do
		end do
	end do
	close (unit4)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_4

!******************************************************************************
       return
       end
