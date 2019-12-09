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
    integer :: ntstep=30
    integer :: ncoef=3
    integer :: nit=100
    !
    integer :: nintf
    ! put equations to be solve here
    type(scalar_t), pointer :: scalar => NULL()
    !type(uvwp_t), pointer :: uvwp => NULL()
    !type(energy_t), pointer :: energy => null()
    ! multiphase solution
    integer :: iphase
    type(phase_t) :: phase(0:NPHASES)
    type(interaction_array_t) :: interaction(NPHASES,NPHASES)
    ! segregated solver coefficients
    real, allocatable, dimension(:) :: ap,anb,b,phic
    ! homogeneous level properties
    type(properties_t) :: prop
  end type
 contains

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
    ! set phase component set
    phys%phase(GAS)%COMPONENT_SET(WAT_VAP)=.true.
    phys%phase(LIQUID)%COMPONENT_SET(WAT_LIQ)=.true.

    call init_properties(phys%phase(GAS)%prop,WAT_VAP,geom%ne)
    call init_properties(phys%phase(LIQUID)%prop,WAT_LIQ,geom%ne)
    call init_properties(phys%phase(MIXTURE)%prop,0,geom%ne)! allocate mixture

! construct equations
    phys%phase(MIXTURE)%uvwp => construct_uvwp(geom,phys%phase(MIXTURE)%prop,phys%dt,MIXTURE)

    do phase=1,NPHASES
      phys%phase(phase)%COMPONENT_SET=.false.

      phys%phase(phase)%uvwp => construct_uvwp(geom,phys%phase(phase)%prop,phys%dt,phase)

      phys%phase(phase)%vfr => construct_vfr(geom,phys%phase(phase)%uvwp%mip,phys%phase(phase)%prop,phase)

      phys%phase(phase)%energy => construct_energy(geom,phys%phase(phase)%uvwp%mip,phys%phase(phase)%prop,phase)
    enddo
    phys%phase(LIQUID)%ndf => construct_ndf(geom,phys%phase(LIQUID)%uvwp%mip,phys%phase(LIQUID)%prop,LIQUID)

! construct source terms for interaction between phases
    call update_mixture_properties(phys%phase,geom%ne)! initialize mixture
    do phase=1,NPHASES
      do phase2=phase+1,NPHASES
        phys%interaction(phase,phase2)%mixture=>phys%phase(0)
        phys%interaction(phase,phase2)%donor_phase=>phys%phase(phase)
        phys%interaction(phase,phase2)%receiver_phase=>phys%phase(phase2)
        phys%interaction(phase2,phase)%mixture=>phys%phase(0)
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
  end subroutine

  subroutine destroy_phys(phys)
    type(phys_t) :: phys

    deallocate(phys%anb,phys%ap,phys%b,phys%phic)

  end subroutine

  subroutine update_time(phys)
    type(phys_t) :: phys
    integer :: P

    !phys%scalar%phi0=phys%scalar%phi
    do p=1,NPHASES
      phys%phase(p)%uvwp%u0=phys%phase(p)%uvwp%u
      phys%phase(p)%uvwp%v0=phys%phase(p)%uvwp%v
      phys%phase(p)%uvwp%w0=phys%phase(p)%uvwp%w
      phys%phase(p)%uvwp%mip0=phys%phase(p)%uvwp%mip

      phys%phase(p)%energy%phi0=phys%phase(p)%energy%phi
      phys%phase(p)%vfr%phi0=phys%phase(p)%vfr%phi
    enddo
    phys%phase(LIQUID)%ndf%phi0=phys%phase(LIQUID)%ndf%phi

  end subroutine

end module
