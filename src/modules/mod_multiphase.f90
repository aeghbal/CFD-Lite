!
!  mod_multiphase.f90
!
! Multiphase module
!
module mod_multiphase
    use mod_properties
    use mod_mfr
    use mod_vfr
    use mod_energy
    implicit none

    type  :: component_t
      ! Store phase and component
      integer :: phase,component
      !
      type(properties_t) :: prop
      !
      type(mfr_t), pointer :: mfr => null()

    end type component_t

    type :: phase_t
      ! Store phase and component
      integer :: phase
      logical :: COMPONENT_SET(NCOMPONENTS)=(/.false.,.false.,.false.,.false.,.false.,.false.,.false./)
      !
      type(properties_t) :: prop
      type(component_t) :: component(NCOMPONENTS)
      !
      type(energy_t), pointer :: energy => null()
      type(vfr_t), pointer :: vfr => null()
      ! DQMOM variables
       real, allocatable, dimension(:) :: a_alpha,b_alpha,Dcoef
       real, allocatable, dimension(:) :: gvolD
       real, allocatable, dimension(:) :: weight,volD

    end type phase_t
  contains

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
            phase_prop(e)=phase_prop(e)+phase%component(c)%mfr%phi(e)*component_prop(e)
          enddo
        end if
      end do
    end subroutine

    subroutine calc_homogeneous_prop(homogeneous_prop,phase,cprop,ne)
      type(phase_t) :: phase(NPHASES)
      type(properties_t) :: homogeneous_properties
      character(len=*) :: cprop
      integer :: p,ne,e
      real, pointer, dimension(:) :: phase_prop,homogeneous_prop

      homogeneous_prop=>homogeneous_properties%get_prop(cprop)
      homogeneous_prop=0.
      do p=1,NPHASES
          phase_prop => phase(p)%prop%get_prop(cprop)
          do e=1,ne
            homogeneous_prop(e)=homogeneous_prop(e)+phase(p)%vfr%phi(e)*phase_prop(e)
          enddo
      end do
    end subroutine

end module mod_multiphase
