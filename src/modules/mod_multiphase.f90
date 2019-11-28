!
!  mod_multiphase.f90
!
! Multiphase module
!
module mod_multiphase
    use mod_properties
    use mod_energy
    use mod_uvwp
    use mod_vfr
    use mod_mfr
    use mod_ndf
    implicit none

    type  :: component_t
      ! Store phase and component
      integer :: phase,component
      !
      type(properties_t) :: prop
      !
      type(mfr_t), pointer :: mfr => null()
      real, allocatable :: vfr(:)

    end type component_t

    type :: phase_t
      ! Store phase and component
      integer :: phase
      logical :: COMPONENT_SET(NCOMPONENTS)
      !
      type(properties_t) :: prop
      type(component_t) :: component(NCOMPONENTS)
      !
      type(uvwp_t), pointer :: uvwp => null()
      type(energy_t), pointer :: energy => null()
      type(vfr_t), pointer :: vfr => null()
      type(ndf_t), pointer :: ndf => null()

      !
      real, allocatable, dimension(:) :: mfr
    end type phase_t

    integer, parameter :: N_INTERACTIONS=4,NUCLEATION=1,MASS_EXCH=2,MOMENTUM_EXCH=3,ENERGY_EXCH=4
    type :: phase_interaction_t
      real :: sgn=1.0
      real, pointer, dimension(:) :: src=>null()
      real, pointer, dimension(:) :: sdc=>null()! deferred correction coefficient of this src
    end type

    type :: interaction_array_t
      type(phase_t) , pointer :: donor_phase=>null()
      type(phase_t) , pointer :: receiver_phase=>null()
      type(phase_interaction_t), dimension(N_INTERACTIONS) :: scheme
     contains
      procedure :: update
    end type
  contains
    subroutine update(interaction,scheme,geom)
      class(interaction_array_t) :: interaction
      type(geometry_t) :: geom
      integer :: scheme

      select case(scheme)
        case(NUCLEATION)
          call calc_nucleation_src(interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(MASS_EXCH)
          call calc_mass_exch_src(interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(MOMENTUM_EXCH)
          call calc_momentum_exch_src(interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(ENERGY_EXCH)
          call calc_energy_exch_src(interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
      end select

    end subroutine

    subroutine calc_nucleation_src(dnr,rcvr,scheme,geom)
      type(phase_interaction_t) :: scheme
      type(phase_t) :: dnr,rcvr
      type(geometry_t) :: geom
      integer :: e
      real :: J,eta,r_crit,DGibbs,m_crit
      real, parameter :: z=1.0,pi=3.14,gam=1.4,L=2230,K_b=1.38e-23,R_u=8.3145,M_wat=18


      do e=1,geom%ne

        DGibbs=dnr%prop%sigma(e)*(dnr%vfr%phi(e)/dnr%ndf%phi(e)*geom%vol(e))
        r_crit=2*dnr%prop%sigma(e)/dnr%prop%rho(e)/DGibbs
        eta=2*(gam-1)/(gam+1)*L/R_u/rcvr%energy%t(e)*(L/R_u/rcvr%energy%t(e)-0.5)
        J=z/(1+eta)*sqrt(2*dnr%prop%sigma(e)/pi/M_wat**3)*rcvr%prop%rho(e)**2/dnr%prop%rho(e)*exp(-4*pi*r_crit**2*dnr%prop%sigma(e)/3/K_b/rcvr%energy%t(e))
        m_crit=4./3.*pi*r_crit**3*rcvr%prop%rho(e)
        scheme%src(e) = m_crit*dnr%vfr%phi(e)*J

      enddo
    end subroutine

    subroutine calc_mass_exch_src(dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr
      type(geometry_t) :: geom

    end subroutine

    subroutine calc_momentum_exch_src(dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr
      type(geometry_t) :: geom
      real :: C_d,beta


    end subroutine

    subroutine calc_energy_exch_src(dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr
      type(geometry_t) :: geom

    end subroutine

    subroutine calc_algebraic_mfr(phase,icomponent,ne)
      type(phase_t) :: phase
      integer :: icomponent,c,ne,e

      phase%component(icomponent)%mfr%phi=1.
      do c=1,NCOMPONENTS
        if(phase%COMPONENT_SET(c) .and. c/=icomponent) then
          do e=1,ne
            phase%component(icomponent)%mfr%phi(e)=phase%component(icomponent)%mfr%phi(e)-phase%component(c)%mfr%phi(e)
          enddo
        end if
      end do

    end subroutine

    subroutine calc_algebraic_vfr(phase,iphase,ne)
      type(phase_t) :: phase(NPHASES)
      integer :: iphase,p,ne,e

      phase(iphase)%vfr%phi=1.
      do p=1,NPHASES
        if(p/=iphase) then
          do e=1,ne
            phase(iphase)%vfr%phi(e)=phase(iphase)%vfr%phi(e)-phase(p)%vfr%phi(e)
          enddo
        end if
      end do

    end subroutine

    subroutine calc_phase_prop(phase,cprop,ne)
      type(phase_t) :: phase
      character(len=*) :: cprop
      integer :: c,ne,e
      real, pointer, dimension(:) :: phase_prop,component_prop

      phase_prop=>phase%prop%get_prop(cprop)
      phase_prop=0.
      do c=1,NCOMPONENTS
        if(phase%COMPONENT_SET(c)) then
          component_prop => phase%component(c)%prop%get_prop(cprop)
          do e=1,ne
            phase_prop(e)=phase_prop(e)+phase%component(c)%vfr(e)*component_prop(e)
          enddo
        end if
      end do
    end subroutine

    subroutine calc_mixture_prop(mixture_prop,phase,cprop,ne)
      type(phase_t) :: phase(NPHASES)
      type(properties_t) :: mixture_properties
      character(len=*) :: cprop
      integer :: p,ne,e
      real, pointer, dimension(:) :: phase_prop,mixture_prop

      mixture_prop=>mixture_properties%get_prop(cprop)
      mixture_prop=0.
      do p=1,NPHASES
          phase_prop => phase(p)%prop%get_prop(cprop)
          do e=1,ne
            mixture_prop(e)=mixture_prop(e)+phase(p)%vfr%phi(e)*phase_prop(e)
          enddo
      end do
    end subroutine

    subroutine calc_phase_field(phase,cprop,ne)
      type(phase_t) :: phase
      integer :: c,ne,e
      real, pointer, dimension(:) :: phase_prop,component_prop
      character(len=*) :: cprop

      phase_prop=>phase%prop%get_prop(cprop)
      phase_prop=0.
      do c=1,NCOMPONENTS
        if(phase%COMPONENT_SET(c)) then
          component_prop => phase%component(c)%prop%get_prop(cprop)
          do e=1,ne
            phase_prop(e)=phase_prop(e)+phase%component(c)%mfr%phi(e)*component_prop(e)
          enddo
        end if
      end do
    end subroutine

    subroutine calc_mixture_field(mixture_eqn,phase,ne)
      type(phase_t) :: phase(NPHASES)
      class(equation_t) :: mixture_eqn
      integer :: p,ne,e
      real, pointer, dimension(:) :: phase_field,mixture_field,phase_u,phase_v,phase_w

      select type(mixture_eqn)
        type is(uvwp_t)
          do e=1,ne
            mixture_eqn%u(e)=0.
            mixture_eqn%v(e)=0.
            mixture_eqn%w(e)=0.
            do p=1,NPHASES
              mixture_eqn%u(e)=mixture_eqn%u(e)+phase(p)%vfr%phi(e)*phase(p)%uvwp%u(e)
              mixture_eqn%v(e)=mixture_eqn%v(e)+phase(p)%vfr%phi(e)*phase(p)%uvwp%v(e)
              mixture_eqn%w(e)=mixture_eqn%w(e)+phase(p)%vfr%phi(e)*phase(p)%uvwp%w(e)
            enddo
          enddo
        class default
          do e=1,ne
            mixture_eqn%phi(e)=0.
            do p=1,NPHASES
              mixture_eqn%phi(e)=mixture_eqn%phi(e)+phase(p)%vfr%phi(e)*phase(p)%uvwp%phi(e)
            enddo
          enddo
      end select

    end subroutine


end module mod_multiphase
