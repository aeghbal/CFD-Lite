
program cfdlite
  use mod_cell
  use mod_physics
  use mod_solver
  use mod_vtu_output
  use mod_vtk
  implicit none

  character(len=180) :: filename,projPath

  integer :: l,i
  type(geometry_t) :: geom
  type(phys_t) :: phys
  integer :: tstep,icoef,it
  real :: residual_rms_i,residual_rms_f,residual_max
  real :: start, finish

#ifdef Catalyst
    print*, "Catalyst Adaptor: Enabled"
#elseif defined VTK
    Print*, "VTK Adaptor: Enabled"
#endif

! this program takes only 1 argument which is the cgns file
  call get_command_argument(1,filename)
  projPath='./'
  l=len_trim(filename)
  do i=l,1,-1
    if(filename(i:i)=='/') exit
  end do
  if(i>0) projPath(1:i)=filename(1:i)

  write(*,*) 'Project path is :'//trim(projPath)
  call cpu_time(start)
! construct the geometry that has only 1 data block and c2b interfaces as many as 2D sections of it
  call cell_input(geom,filename,phys%n_subdomains)
!
! construct physics
  call construct_physics(phys,geom)

#ifdef Catalyst
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cfd_catalyst( 1, phys%uvwp%p(1), phys%uvwp%gpc(1), phys%uvwp%u(1), phys%uvwp%v(1), phys%uvwp%w(1), phys%energy%t(1), geom%x(1), geom%y(1), geom%z(1), geom%mg%e2vx(1),&
						   geom%nvx,geom%ne, geom%mg%esec(1,1), geom%mg%etype(1),geom%mg%nsec,geom%mg%ne2vx_max, 0, &
						   phys%ntstep, phys%dt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

!
  write(*,'(5x,A16,x,A5,x,A15,A9,3x,A9,3x,A9)') 'solve eqn       ','nit','residual(rms)','initial','final','max'
  do tstep=1,phys%ntstep
    do icoef=1,phys%ncoef
      call update_boundaries(phys,geom)
      ! solve scalar
      !do phys%iphase=1,NPHASES
        !call solve_scalar(phys%scalar,phys%prop,geom,phys%dt,phys%nit,phys%ap,phys%anb,phys%b,phys%phic)
        call solve_uvwp(phys%uvwp,phys%prop,geom,phys%dt,phys%nit,phys%ap,phys%anb,phys%b,phys%phic,phys%subdomain,phys%intf,phys%n_subdomains)

        !call solve_energy(phys%energy,phys%prop,geom,phys%dt,phys%nit,phys%ap,phys%anb,phys%b,phys%phic)
      !enddo
    end do

    ! update time level
    call update_time(phys)
    write(*,'(A,i5,A)') '------------------------------------time step(',tstep,')'
#ifdef Catalyst
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call cfd_catalyst( 2, phys%uvwp%p(1), phys%uvwp%gpc(1), phys%uvwp%u(1), phys%uvwp%v(1), phys%uvwp%w(1), phys%energy%t(1), geom%x(1), geom%y(1), geom%z(1), geom%mg%e2vx(1),&
						   geom%nvx,geom%ne, geom%mg%esec(1,1), geom%mg%etype(1),geom%mg%nsec,geom%mg%ne2vx_max, tstep, &
						   phys%ntstep, phys%dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

#ifdef VTK
	if(mod(tstep,10)==0) then
    	call vtk_data(phys%uvwp,phys%energy,geom,tstep)
	endif
#else
    if(mod(tstep,10)==0) then
      call write_vtubin(phys%uvwp,geom,projPath,tstep)
     ! call write_vtubin(phys%energy,geom,projPath,tstep)
    endif
#endif



  end do

  !write_vtu
  call write_vtubin(phys%uvwp,geom,projPath,tstep)

#ifdef Catalyst
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cfd_catalyst( 3, phys%uvwp%p(1), phys%uvwp%gpc(1), phys%uvwp%u(1), phys%uvwp%v(1), phys%uvwp%w(1), phys%energy%t(1), geom%x(1), geom%y(1), geom%z(1), geom%mg%e2vx(1),&
						   geom%nvx,geom%ne, geom%mg%esec(1,1), geom%mg%etype(1),geom%mg%nsec,geom%mg%ne2vx_max, tstep, &
						   phys%ntstep, phys%dt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

 call cpu_time(finish)
    print '("Time (total) = ",f10.5," seconds; Time (Pressure Correction Eqn) = ",f10.5," seconds.")',finish-start, geom%elapt
  !call destroy_scalar(phys%scalar)
  call destroy_phys(phys)
  call destroy_geom(geom)

end program

