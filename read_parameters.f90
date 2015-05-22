subroutine read_parameters

! This subroutine can read the data on the mesh and flow conditions.

       use common_block

       implicit none
       
	character (len = 64) :: output_file_name_1, str, word, tmp
	integer :: unit1

	output_file_name_1 = 'cube.parameters.dat'
	unit1 = 1
	open (unit=unit1,file=output_file_name_1)
	read (unit1,'(a)') !'Parameters:'
	read (unit1, '(a)') !' x: '
	read (unit1, '(g14.6)') LX
	read (unit1, '(a)') !' y: '
	read (unit1, '(g14.6)') LY
	read (unit1, '(a)') !' z: '
	read (unit1, '(g14.6)') LZ
	read (unit1, '(a)') !' gapx: '
	read (unit1, '(g14.6)') gap
	read (unit1, '(a)') !' gapy: '
	read (unit1, '(g14.6)') !gapy
	read (unit1, '(a)') !' gapz: '
	read (unit1, '(g14.6)') !gapz
	read (unit1, '(a)') !' xtrans: '
	read (unit1, '(g14.6)') xtrans
	read (unit1, '(a)') !' ytrans: '
	read (unit1, '(g14.6)') ytrans
	read (unit1, '(a)') !' ztrans: '
	read (unit1, '(g14.6)') ztrans

	close (unit1)
	write (*, '(a)') ' '
	write (*, '(a)' , advance = "no") ' File opened: '
	write (*, '(a)') output_file_name_1
	       
       return
end
