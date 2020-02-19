module mod_physics
  use mod_cell
  use mod_util
  use mod_eqn_setup
  use mod_scalar
  use mod_energy
  use mod_multiphase
  use mod_uvwp
  use mod_subdomains
  implicit none
  ! Physics

  type :: phys_t
    ! temporal properties
    real :: dt=0.01
    integer :: ntstep=1000
    integer :: ncoef=3
    integer :: nit=100
    integer :: n_subdomains=4
    !
    integer :: nintf_c2b
    ! put equations to be solve here
    type(scalar_t), pointer :: scalar => NULL()
    type(uvwp_t), pointer :: uvwp => NULL()
    type(energy_t), pointer :: energy => null()
    ! multiphase solution
    integer :: iphase
    type (phase_t) :: phase(NPHASES)
    ! segregated solver coefficients
    real, allocatable, dimension(:) :: ap,anb,b,phic
    type(subdomain_t), allocatable :: subdomain(:)
    type(intf_t) , allocatable :: intf(:,:)
    ! homogeneous level properties
    type(properties_t) :: prop
  end type
 contains

   subroutine update_boundaries(phys,geom)
    implicit none
    type(phys_t) :: phys
    type(geometry_t) :: geom
    integer :: i

    do i=1,phys%nintf_c2b
      !call phys%scalar%bcs(i)%coef(geom,phys%scalar)
      call phys%uvwp%bcs(i)%coef(geom,phys%uvwp,phys%prop)
      call phys%energy%bcs(i)%coef(geom,phys%energy,phys%prop)
    end do

  end subroutine

  subroutine construct_physics(phys,geom)
    use mod_solver_cuda
    implicit none
    type(phys_t) :: phys
    type(geometry_t) :: geom
    integer :: length,c

    allocate(phys%ap(geom%ne));phys%ap=0.
    allocate(phys%b(geom%ne));phys%b=0.
    length = 2*geom%nf-geom%nbf
    allocate(phys%anb(length));phys%anb=0.
    length = geom%ne+geom%nbf
    allocate(phys%phic(length));phys%phic=0.
    phys%nintf_c2b=geom%mg%nintf_c2b
! construct subdomains
    if(phys%n_subdomains>1) then
      call construct_subdomains(phys%subdomain,phys%n_subdomains,phys%intf,geom)
!$ifdef CUDA
      do c=1,phys%n_subdomains
        call allocate_device_data(phys%subdomain(c))
      enddo
!$endif
    endif

! initialize properties
    call init_properties(phys%prop,0,geom%ne)
! construct equations
    phys%uvwp => construct_uvwp(geom,phys%prop,phys%dt)

    phys%energy => construct_energy(geom,phys%uvwp%mip,phys%prop)

  end subroutine

  subroutine destroy_phys(phys)
    type(phys_t) :: phys
    integer :: c,cnb

    deallocate(phys%anb,phys%ap,phys%b,phys%phic)
    if(phys%n_subdomains>1) then
      do c=1,phys%n_subdomains
        deallocate(phys%subdomain(c)%ap)
        deallocate(phys%subdomain(c)%b)
        deallocate(phys%subdomain(c)%anb)
        deallocate(phys%subdomain(c)%phic)
        deallocate(phys%subdomain(c)%ef2nb)
        deallocate(phys%subdomain(c)%ef2nb_idx)
        do cnb=1,phys%n_subdomains
          deallocate(phys%intf(c,cnb)%index1)
        end do
      end do
      deallocate(phys%subdomain)
      deallocate(Phys%intf)
    end if


  end subroutine

  subroutine update_time(phys)
    type(phys_t) :: phys

    !phys%scalar%phi0=phys%scalar%phi

    phys%uvwp%u0=phys%uvwp%u
    phys%uvwp%v0=phys%uvwp%v
    phys%uvwp%w0=phys%uvwp%w
    phys%uvwp%mip0=phys%uvwp%mip

    phys%energy%phi0=phys%energy%phi
  end subroutine

end module
