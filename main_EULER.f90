PROGRAM main_EULER

!                           Euler - 3D Solver
!
!
!      THIS PROGRAM SOLVES THE EULER EQUATIONS FOR THREE DIMENSIONS 
!      USING AN THETRAHEDRAL MESH WITH CELL CORNER STORAGE OF THE VARIABLES.
!

!         The code was developed by Tatiele D. Ferreira PhD' student.
!                       University of Campinas
!                  Faculty of Chemical Engineering
!
!                 Campinas, Sao Paulo August, 2014.

       use common_block

!  Open the file "euler.log" which is used to plot the convergence
!  of the calculation via the separate program "PLTCONV".

      OPEN(UNIT=3,FILE='euler.log')
!
!  Call subroutine "read_data" to read in the flow conditions.
!
      call read_data
!
      call read_parameters
!
       print*, lx, ly, lz, gap
!
       print*, xtrans, ytrans, ztrans

       DNI=nint(lx/gap)
       DNJ=nint(ly/gap)
       DNK=nint(lz/gap)

       NI = DNI + 1.
       NJ = DNJ + 1.
       NK = DNK + 1.

       NPOINTS=NI*NJ*NK
!
      call read_porosity
!
      call check_porosity
!
!    Call subroutine "mesh" to set up the grid coordinates,
!    element areas and projected lengths of the sides of the elements.
!
!     Call subroutine "check_grid" to check that the areas and projected
!     lengths are correct.
!
      call check_grid
!
!
       allocate ( ROVX(1:NI,1:NJ,1:NK) )
       allocate ( ROVY(1:NI,1:NJ,1:NK) )
       allocate ( ROVZ(1:NI,1:NJ,1:NK) )
       allocate ( ROE(1:NI,1:NJ,1:NK) )
       allocate ( RO(1:NI,1:NJ,1:NK) )

       allocate ( P(1:NI,1:NJ,1:NK) )
       allocate ( HSTAG(1:NI,1:NJ,1:NK) )
!
!    Call subroutine "inguess" this gives some arbitrary values to the
!    flow variables so that you can plot the grid. These values are not used
!    at all in the subsequent calculations.
!
      call inguess	
!
!     It is possible to call subroutine "OUTPUT" here to plot out the grid
!     you have generated.
!
!     Subroutine  "crude_guess" is what its name says. It enables  to
!     start a calculation and obtain a solution but it will take longer than
!     necessary. When the program is working one should replace it
!     with "flow_guess" to obtain a better guess and a faster solution.
!
       allocate ( VX(1:NI,1:NJ,1:NK) )
       allocate ( VY(1:NI,1:NJ,1:NK) )
       allocate ( VZ(1:NI,1:NJ,1:NK) )

       allocate ( ROVX_OLD(1:NI,1:NJ,1:NK) )
       allocate ( ROVY_OLD(1:NI,1:NJ,1:NK) )
       allocate ( ROVZ_OLD(1:NI,1:NJ,1:NK) )
       allocate ( ROE_OLD(1:NI,1:NJ,1:NK) )
       allocate ( RO_OLD(1:NI,1:NJ,1:NK) )
       allocate ( diffrovx(1:NI,1:NJ,1:NK) )
!
        call crude_guess
!
!        call flow_guess
!
!      It is possible to  call "output" here to plot out your initial guess of
!      the flow field.
!
!     Call subroutine "set_timestep" to set the length of the timestep.
!     Initially this is a constant time step based on a conservative guess
!     of the Mach number.
!
      call set_timestep
!
       allocate ( FLUXI_MA(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( FLUXJ_MA(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( FLUXK_MA(1:(NI-1),1:(NJ-1),1:NK) )
       allocate ( FLOW(1:NI) )

       allocate ( FLUXI_XMOM(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( FLUXJ_XMOM(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( FLUXK_XMOM(1:(NI-1),1:(NJ-1),1:NK) )

       allocate ( FLUXI_YMOM(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( FLUXJ_YMOM(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( FLUXK_YMOM(1:(NI-1),1:(NJ-1),1:NK) )

       allocate ( FLUXI_ZMOM(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( FLUXJ_ZMOM(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( FLUXK_ZMOM(1:(NI-1),1:(NJ-1),1:NK) )

       allocate ( FLUXI_ENTH(1:NI,1:(NJ-1),1:(NK-1)) )
       allocate ( FLUXJ_ENTH(1:(NI-1),1:NJ,1:(NK-1)) )
       allocate ( FLUXK_ENTH(1:(NI-1),1:(NJ-1),1:NK) )

       allocate ( DELRO(1:NI,1:NJ,1:NK) )
       allocate ( DELROVX(1:NI,1:NJ,1:NK) )
       allocate ( DELROVY(1:NI,1:NJ,1:NK) )
       allocate ( DELROVZ(1:NI,1:NJ,1:NK) )
       allocate ( DELROE(1:NI,1:NJ,1:NK) )

       allocate ( ROINLET(1:NJ,1:NK) )
       allocate ( T(1:NI,1:NJ,1:NK) )
       allocate ( V(1:NI,1:NJ,1:NK) )
!
       print*,'START THE ITERATION PROCESS'
!
!
!************************************************************************
!      START THE TIME STEPPING DO LOOP FOR "NSTEPS" LOOPS.
!************************************************************************
!
 !     open(unit=200,file='FLUX.dat')
 !     write(200,"(A)") 'NSTEP fluxi_mass fluxJ_mass fluxK_mass'
!
      DO 1000 NSTEP = 1,NSTEPS
!
        call set_others

!
!       print*,'here9'
!
      call inout_bc
!
!
      call flux_ma
!
      call wall_bc
!
      call flux_energy
!

      call sum_fluxes (fluxi_ma,fluxj_ma,fluxk_ma,ro,delro)
      call sum_fluxes (fluxi_enth,fluxj_enth,fluxk_enth,roe,delroe)
!

      call smooth(ro)
      call smooth(roe)

! --------------------------------------------------------------------------

      call flux_mom
!
      call sum_fluxes (fluxi_xmom,fluxj_xmom,fluxk_xmom,rovx,delrovx)
      call sum_fluxes (fluxi_ymom,fluxj_ymom,fluxk_ymom,rovy,delrovy)
      call sum_fluxes (fluxi_zmom,fluxj_zmom,fluxk_zmom,rovz,delrovz)
!
      call smooth(rovx)
      call smooth(rovy)
      call smooth(rovz)
!
!-------------------------------------------------------------------------
!
!       print*,'here13'
!
!     CHECK CONVERGENCE AND WRITE OUT SUMMARY EVERY 5 STEPS
!
      IF(MOD(NSTEP,5).EQ.0) THEN
!
	 call check_conv
 !
      ENDIF
!
!
      IF(MOD(NSTEP,NDUMP).EQ.0) THEN 
            CALL output(NSTEP)
      ENDIF
!
!     GO TO 2000 IF CONVERGED TO THE INPUT TOLERANCE  "conlim"
!
      IF(EMAX.LT.CONLIM.AND.EAVG.LT.(0.5*CONLIM)) GO TO 2000
!
!     END OF THE MAIN TIME STEPPING DO LOOP
!
 1000 CONTINUE
!
!******************************************************************************
!
 2000 CONTINUE
!
      IF(EMAX.LT.CONLIM.AND.EAVG.LT.(0.5*CONLIM)) THEN
	  WRITE(6,*) ' CALCULATION CONVERGED IN ',NSTEP,' ITERATIONS'
	  WRITE(6,*) ' TO A CONVERGENCE LIMIT OF ', CONLIM
      ENDIF
!
!   CALCULATION FINISHED. CALL "output" TO WRITE THE PLOTTING FILE.
!
      call output(NSTEP)
!
!      close(200)

       stop
       END
