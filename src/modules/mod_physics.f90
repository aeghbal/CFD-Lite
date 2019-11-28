module mod_physics
  use mod_cell
  use mod_util
  use mod_eqn_setup
  use mod_scalar
  use mod_energy
  use mod_multiphase
  use mod_uvwp
  use mod_vfr
  use mod_mfr
  use mod_ndf
  implicit none
  ! Physics
  type :: phys_t
    ! temporal properties
    real :: dt=0.01
    integer :: ntstep=100
    integer :: ncoef=3
    integer :: nit=100
    !
    integer :: nintf
    ! put equations to be solve here
    type(scalar_t), pointer :: scalar => NULL()
    type(uvwp_t), pointer :: uvwp => NULL()
    type(energy_t), pointer :: energy => null()
    ! multiphase solution
    integer :: iphase
    type(phase_t) :: phase(NPHASES)
    type(interaction_array_t) :: interaction(NPHASES,NPHASES)
    ! segregated solver coefficients
    real, allocatable, dimension(:) :: ap,anb,b,phic
    ! homogeneous level properties
    type(properties_t) :: prop
  end type
 contains

   subroutine update_boundaries(phys,geom)
    implicit none
    type(phys_t) :: phys
    type(geometry_t) :: geom
    integer :: i

    do i=1,phys%nintf
      !call phys%scalar%bcs(i)%coef(geom,phys%scalar)
      call phys%uvwp%bcs(i)%coef(geom,phys%uvwp,phys%prop)
      call phys%energy%bcs(i)%coef(geom,phys%energy,phys%prop)
    end do

  end subroutine

  subroutine construct_physics(phys,geom)
    implicit none
    type(phys_t), target :: phys
    type(geometry_t) :: geom
    integer :: length,phase,phase2,i

    allocate(phys%ap(geom%ne));phys%ap=0.
    allocate(phys%b(geom%ne));phys%b=0.
    length = 2*geom%nf-geom%nbf
    allocate(phys%anb(length));phys%anb=0.
    length = geom%ne+geom%nbf
    allocate(phys%phic(length));phys%phic=0.
    phys%nintf=geom%mg%nintf

! initialize properties
    call init_properties(phys%prop,0,geom%ne)
! construct equations
    do phase=1,NPHASES
      phys%phase(phase)%COMPONENT_SET=.false.

      phys%phase(phase)%vfr => construct_vfr(geom,phys%phase(phase)%uvwp%mip,phys%prop)

      phys%phase(phase)%uvwp => construct_uvwp(geom,phys%prop,phys%dt,phys%phase(phase)%vfr%phi,phys%phase(phase)%vfr%phi0)

      phys%phase(phase)%energy => construct_energy(geom,phys%phase(phase)%uvwp%mip,phys%prop,phys%phase(phase)%vfr%phi,phys%phase(phase)%vfr%phi0)
    enddo

    phys%phase(LIQUID)%ndf => construct_ndf(geom,phys%phase(phase)%uvwp%mip,phys%prop,phys%phase(phase)%vfr%phi,phys%phase(phase)%vfr%phi0)
    ! construct source terms for interaction between phases
    do phase=1,NPHASES
      do phase2=phase+1,NPHASES
        phys%interaction(phase,phase2)%donor_phase=>phys%phase(phase)
        phys%interaction(phase,phase2)%receiver_phase=>phys%phase(phase2)
        phys%interaction(phase2,phase)%donor_phase=>phys%phase(phase2)
        phys%interaction(phase2,phase)%receiver_phase=>phys%phase(phase)
        do i=1,N_INTERACTIONS
          ! receiver-donor
          phys%interaction(phase,phase2)%scheme(i)%sgn=1.0
          allocate(phys%interaction(phase,phase2)%scheme(i)%src(geom%ne));phys%interaction(phase,phase2)%scheme(i)%src=0.
          allocate(phys%interaction(phase,phase2)%scheme(i)%sdc(geom%ne));phys%interaction(phase,phase2)%scheme(i)%sdc=0.
          ! donor-receiver points to above arrays with negative sign
          phys%interaction(phase2,phase)%scheme(i)%sgn=-1.0
          phys%interaction(phase2,phase)%scheme(i)%src=>phys%interaction(phase,phase2)%scheme(i)%src
          phys%interaction(phase2,phase)%scheme(i)%sdc=>phys%interaction(phase,phase2)%scheme(i)%sdc
        end do
      end do
    end do
    ! set phase component set
    phys%phase(GAS)%COMPONENT_SET(WAT_VAP)=.true.
    phys%phase(LIQUID)%COMPONENT_SET(WAT_LIQ)=.true.
!
  end subroutine

  subroutine destroy_phys(phys)
    type(phys_t) :: phys

    deallocate(phys%anb,phys%ap,phys%b,phys%phic)

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
