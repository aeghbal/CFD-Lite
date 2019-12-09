
program cfdite
  use mod_cell
  use mod_task_launcher
  use mod_vtu_output
  implicit none

  character(len=180) :: filename,projPath

  integer :: l,i,p
  type(geometry_t) :: geom
  type(phys_t) :: phys
  integer :: tstep,icoef,it
  real :: residual_rms_i,residual_rms_f,residual_max

! this program takes only 1 argument which is the cgns file
  call get_command_argument(1,filename)
  projPath='./'
  l=len_trim(filename)
  do i=l,1,-1
    if(filename(i:i)=='/') exit
  end do
  if(i>0) projPath(1:i)=filename(1:i)

  write(*,*) 'Project path is :'//trim(projPath)
! construct the geometry that has only 1 data block and c2b interfaces as many as 2D sections of it
  call cell_input(geom,filename)
!
! construct physics
  call construct_physics(phys,geom)
!
  write(*,'(5x,A16,x,A5,x,A15,A9,3x,A9,3x,A9)') 'solve eqn       ','nit','residual(rms)','initial','final','max'
  do tstep=1,phys%ntstep
    do icoef=1,phys%ncoef

      ! solve continuity/mass
      call solve_vfr(phys,geom)
      call solve_properties(phys,geom)
      ! solve momentum
      call solve_uvw(phys,geom)
      !
      call solve_pressure(phys,geom)
      ! update properties
      !call update_properties(phys,geom)

      if(icoef>1) then
        ! solve energy
        call solve_energy(phys,geom)
        ! solve population balance eqn
        call solve_ndf(phys,geom)
        ! solve mass fraction
        !call solve_mfr(phys,geom)
        ! update properties
        call solve_properties(phys,geom)

      endif
    end do

    ! update time level
    call update_time(phys)
    write(*,'(A,i5,A)') '------------------------------------time step(',tstep,')'

    if(mod(tstep,10)==0) then
      do p=0,NPHASES
        if(p>0) call write_vtubin(phys%phase(p)%vfr,geom,phys%phase(p)%prop,projPath,tstep)
        call write_vtubin(phys%phase(p)%uvwp,geom,phys%phase(p)%prop,projPath,tstep)
        if(p>0)call write_vtubin(phys%phase(p)%energy,geom,phys%phase(p)%prop,projPath,tstep)
      enddo
      call write_vtubin(phys%phase(LIQUID)%ndf,geom,phys%phase(LIQUID)%prop,projPath,tstep)
    endif

  end do

  !write_vtu
  !call write_vtubin(phys%scalar,geom,projPath,1)

  !call destroy_scalar(phys%scalar)
  call destroy_phys(phys)
  call destroy_geom(geom)

end program

