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

    end type phase_t

    integer, parameter :: N_INTERACTIONS=4,NUCLEATION=1,MASS_EXCH=2,MOMENTUM_EXCH=3,ENERGY_EXCH=4
    type :: phase_interaction_t
      real :: sgn=1.0
      real, pointer, dimension(:) :: src=>null()
      real, pointer, dimension(:) :: sdc=>null()! deferred correction coefficient of this src
    end type

    type :: interaction_array_t
      type(phase_t) , pointer :: mixture=>null()
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
return
      select case(scheme)
        case(NUCLEATION)
          call calc_nucleation_src(interaction%mixture,interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(MASS_EXCH)
          call calc_mass_exch_src(interaction%mixture,interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(MOMENTUM_EXCH)
          call calc_momentum_exch_src(interaction%mixture,interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
        case(ENERGY_EXCH)
          call calc_energy_exch_src(interaction%mixture,interaction%donor_phase,interaction%receiver_phase,interaction%scheme(scheme),geom)
      end select

    end subroutine

    subroutine calc_nucleation_src(mixtr,dnr,rcvr,scheme,geom)
      type(phase_interaction_t) :: scheme
      type(phase_t) :: dnr,rcvr,mixtr
      type(geometry_t) :: geom
      integer :: e
      real :: J,eta,r_crit,DGibbs,m_crit,v_m,p_sat
      real, parameter :: z=1.0,pi=3.14,gam=1.4,L=2230,K_b=1.38e-23,R_u=8.3145,M_wat=0.018

! receiver is liquid water
! donor is water vapour
      do e=1,geom%ne
        v_m=M_wat/rcvr%prop%rho(e)
        ! Antoine eqn
        p_sat=10.**(8.1-1770./(240.+dnr%energy%t(e)))
        r_crit=2*rcvr%prop%sigma(e)*v_m/(R_u*dnr%energy%t(e)*log(mixtr%uvwp%p(e)/p_sat))
        eta=2*(gam-1)/(gam+1)*L/R_u/dnr%energy%t(e)*(L/R_u/dnr%energy%t(e)-0.5)
        DGibbs=4./3.*pi*r_crit**2*rcvr%prop%sigma(e)
        J=z/(1+eta)*sqrt(2*rcvr%prop%sigma(e)/pi/M_wat**3)*dnr%prop%rho(e)**2/rcvr%prop%rho(e)*exp(-DGibbs/K_b/dnr%energy%t(e))
        m_crit=4./3.*pi*r_crit**3*rcvr%prop%rho(e)
        scheme%src(e) = m_crit*dnr%vfr%phi(e)*J

      enddo
    end subroutine

    subroutine calc_mass_exch_src(mixtr,dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr,mixtr
      type(geometry_t) :: geom
      integer :: e
      real :: qc,r_d,Nu_c,kn
      real, parameter :: z=1.0,pi=3.14,gam=1.4,L=2230,K_b=1.38e-23,R_u=8.3145,M_wat=0.018,d_wat=2.75e-10
! receiver is liquid water
! donor is water vapour
      ! small droplet model i.e. droplet temp are homogeneous in the droplet
      do e=1,geom%ne
        r_d=(3.*dnr%vfr%phi(e)/4./pi/max(1.,dnr%ndf%phi(e)))**(1./3.)! average droplet radius
        Kn=K_b*dnr%energy%t(e)/1.4/pi/d_wat**2/mixtr%uvwp%p(e)/r_d! Knudsen number based on Boltzmann Gas model
        Nu_c=2./(1.+3.18*Kn)! vapor Nusselt number by Gyarmathy
        qc=dnr%prop%tc(e)/2./r_d*Nu_c*(rcvr%energy%t(e)-dnr%energy%t(e))
        scheme%src(e) = qc/(dnr%energy%phi(e)-rcvr%energy%phi(e))
      end do

    end subroutine

    subroutine calc_momentum_exch_src(mixtr,dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr,mixtr
      type(geometry_t) :: geom
      real :: C_d,beta,Re_d,vel,r_d
      real, parameter :: z=1.0,pi=3.14,gam=1.4,L=2230,K_b=1.38e-23,R_u=8.3145,M_wat=0.018,d_wat=2.75e-10
      integer :: e
! receiver is liquid water
! donor is water vapour
      do e=1,geom%ne
        r_d=(3.*dnr%vfr%phi(e)/4./pi/dnr%ndf%phi(e))**(1./3.)! average droplet radius
        vel=sqrt(rcvr%uvwp%u(e)**2+rcvr%uvwp%v(e)**2+rcvr%uvwp%w(e)**2)
        Re_d=rcvr%prop%rho(e)*vel*2*r_d/rcvr%prop%mu(e)
        if(Re_d<1000.) then
          C_d=24./Re_d*(1.+0.15*Re_d**0.687)! Schiller-Naumann Drag for 0.1<Re<1000
        else
          C_d=0.44
        endif
        beta=3.*dnr%vfr%phi(e)/r_d! droplet surface area density
        scheme%src(e) = 0.
        scheme%sdc(e) = -C_d/8.*beta*dnr%prop%rho(e)

      end do


    end subroutine

    subroutine calc_energy_exch_src(mixtr,dnr,rcvr,scheme,geom)
      type(phase_interaction_t):: scheme
      type(phase_t) :: dnr,rcvr,mixtr
      type(geometry_t) :: geom
      integer :: e

      do e=1,geom%ne
        scheme%src(e) = 0.
      enddo
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
      type(phase_t) :: phase(0:NPHASES)
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

    subroutine calc_phase_mfr(phase,geom)
      type(phase_t) :: phase(0:NPHASES)
      type(geometry_t) :: geom
      integer :: p,e
      real :: mfr0

      do e=1,geom%ne
        mfr0=0.
        do p=1,NPHASES
          mfr0=mfr0+phase(p)%prop%rho(e)*phase(p)%vfr%phi(e)
        enddo
        do p=1,NPHASES
          phase(p)%vfr%mfr(e)= phase(p)%prop%rho(e)*phase(p)%vfr%phi(e)/mfr0
        end do
      end do
    end subroutine

    subroutine calc_component_vfr(vfr,mfr,geom,comp_prop,phase_prop)
      real :: mfr(*),vfr(*)
      type(geometry_t) :: geom
      type(properties_t) :: comp_prop,phase_prop
      integer :: p,e

      do e=1,geom%ne
        vfr(e)= phase_prop%rho(e)*mfr(e)/comp_prop%rho(e)
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
            phase_prop(e)=phase_prop(e)+phase%component(c)%mfr%vfr(e)*component_prop(e)
          enddo
        end if
      end do
    end subroutine

    subroutine calc_mixture_prop(phase,cprop,ne)
      type(phase_t) :: phase(0:NPHASES)
      character(len=*) :: cprop
      integer :: p,ne,e
      real, pointer, dimension(:) :: phase_prop,mixture_prop

      mixture_prop=>phase(MIXTURE)%prop%get_prop(cprop)
      mixture_prop=0.
      do p=1,NPHASES
          phase_prop => phase(p)%prop%get_prop(cprop)
          do e=1,ne
            mixture_prop(e)=mixture_prop(e)+phase(p)%vfr%phi(e)*phase_prop(e)
          enddo
      end do
    end subroutine

    subroutine update_mixture_properties(phase,ne)
      type(phase_t) :: phase(0:NPHASES)
      integer :: p,ne,e

        call calc_mixture_prop(phase,'rho',ne)
        call calc_mixture_prop(phase,'mu',ne)
        call calc_mixture_prop(phase,'cp',ne)
        call calc_mixture_prop(phase,'tc',ne)

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

    subroutine calc_mixture_field(mixture_eqn,phase,geom)
      type(phase_t) :: phase(0:NPHASES)
      type(geometry_t) :: geom
      class(equation_t) :: mixture_eqn
      integer :: p,e,f,enb,lf,idx,lfnb
      real :: r

      select type(mixture_eqn)
        type is(uvwp_t)
          do e=1,geom%ne
            mixture_eqn%u(e)=0.
            mixture_eqn%v(e)=0.
            mixture_eqn%w(e)=0.
            mixture_eqn%dc(e)=0.
            do p=1,NPHASES
              mixture_eqn%u(e)=mixture_eqn%u(e)+phase(p)%vfr%mfr(e)*phase(p)%uvwp%u(e)
              mixture_eqn%v(e)=mixture_eqn%v(e)+phase(p)%vfr%mfr(e)*phase(p)%uvwp%v(e)
              mixture_eqn%w(e)=mixture_eqn%w(e)+phase(p)%vfr%mfr(e)*phase(p)%uvwp%w(e)
              mixture_eqn%dc(e)=mixture_eqn%dc(e)+phase(p)%vfr%mfr(e)*phase(p)%uvwp%dc(e)
            enddo
          enddo

          do f=1,geom%nf
            mixture_eqn%mip(f)=0.
            do p=1,NPHASES
              mixture_eqn%mip(f)=mixture_eqn%mip(f)+phase(p)%vfr%mfr(e)*phase(GAS)%uvwp%mip(f)
            enddo
          end do
        class default
          do e=1,geom%ne
            mixture_eqn%phi(e)=0.
            do p=1,NPHASES
              mixture_eqn%phi(e)=mixture_eqn%phi(e)+phase(p)%vfr%mfr(e)*phase(p)%uvwp%phi(e)
            enddo
          enddo
      end select

    end subroutine

end module mod_multiphase
