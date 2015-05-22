subroutine read_data

! This subroutine can read the data on the mesh and flow conditions.

       use common_block

       implicit none
!
      real   SMOOTH_FAC_IN,CONLIM_IN
!
!  Assign Fortran unit 1 to file 'flow'

       OPEN(UNIT=1,FILE='flow')
!
!NPOINTS represents the number of points
!
!     Read in the flow data from unit 1.
!     You should read in the following variables sequentially
!     RGAS,GAMMA
!     PSTAGIN,TSTAGIN,ALPHA1,PDOWN
!     CFL,SMOOTH_FAC_IN
!     NSTEPS,CONLIM_IN
!
       READ(1,*)
       READ(1,*) RGAS,GAMMA,CD
       READ(1,*)
       READ(1,*) PSTAGIN,TSTAGIN,PDOWN
       READ(1,*)
       READ(1,*) CFL,SMOOTH_FAC_IN
       READ(1,*) 
       READ(1,*) NSTEPS,CONLIM_IN,NDUMP
!
       print*, CD


      open(unit=200,file='READ_DATA.dat')
             write(200,"(A)") 'RGAS GAMMA PSTAGIN TSTAGIN PDOWN CFL SFIN NSTEPS CONIN  NDUMP'
             write(200,"(7F12.5,I10,F12.5,I10)") RGAS,GAMMA,PSTAGIN,TSTAGIN,PDOWN,CFL,&
                                              SMOOTH_FAC_IN,NSTEPS,CONLIM_IN,NDUMP

      close(200)

!
!     Set some other variables that will be used throughout the
!     calculation. Scale the smoothing factor
!     and convergence limits by CFL.

       EMAX       = 1000000.
       EAVG       = EMAX
       CP         = RGAS*GAMMA/(GAMMA-1.)
       CV         = CP/GAMMA
       FGA        = (GAMMA - 1.)/GAMMA
       SMOOTH_FAC = SMOOTH_FAC_IN*CFL
       CONLIM     = CONLIM_IN*CFL


       CLOSE(1)
!
       return
       end
