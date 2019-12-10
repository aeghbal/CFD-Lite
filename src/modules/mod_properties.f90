!
module mod_properties
    real, parameter :: universal_gas_constant=8.314

    integer, parameter :: NPHASES=2,NCOMPONENTS=2
    integer, parameter :: WAT_VAP=1,WAT_LIQ=2,AIR=3
    integer, parameter :: MIXTURE=0,GAS=1,LIQUID=2

    type  :: properties_t
      ! Store phase and component
      integer :: phase,component

      ! Equation of state (ideal gas)
      real :: rgas  ! gas constant
      real :: mw    ! molecular weight

      ! Saturation pressure (Antoine equation)
      real :: pvap_coef(3)
      real, allocatable, dimension(:) :: pvap  ! Vapor pressure
      ! Density
      real :: rho_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: rho
      real, allocatable, dimension(:) :: rho0

      ! Viscosity
      real :: mu_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: mu
      real, allocatable, dimension(:) :: mu0

      ! Specific Heat
      real :: cp_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: cp
      real, allocatable, dimension(:) :: cp0

      ! Thermal Conductivity
      real :: tc_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: tc
      real, allocatable, dimension(:) :: tc0

      ! Diffusion Coeficient
      real :: diff_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: diff
      real, allocatable, dimension(:) :: diff0

      ! Surface Tension
      real :: sigma_coef(5) ! Polynomial coeficients
      real, allocatable, dimension(:) :: sigma
      real, allocatable, dimension(:) :: sigma0
     contains
      procedure :: get_prop
    end type

  contains

    function get_prop(property,cprop) result(prop)
      class(properties_t), target :: property
      real, pointer, dimension(:) :: prop
      character(len=*) :: cprop

      select case(trim(cprop))
        case('rho')
          prop=>property%rho
        case('mu')
          prop=>property%mu
        case('cp')
          prop=>property%cp
        case('tc')
          prop=>property%tc
        case('diff')
          prop=>property%diff
        case('sigma')
          prop=>property%sigma
      end select
    end function

  subroutine init_properties(prop,icomponent,ne)
    implicit none
    type (properties_t) :: prop
    integer :: icomponent,ne
    !
    allocate(prop%rho(ne))
    allocate(prop%rho0(ne))

    allocate(prop%mu(ne))
    allocate(prop%cp(ne))
    allocate(prop%tc(ne))
    allocate(prop%sigma(ne))
    allocate(prop%pvap(ne))
    select case(icomponent)
      case(WAT_LIQ)
        prop%rho=997
        prop%rho0=prop%rho
        prop%mu=0.001
        prop%cp=4200
        prop%tc=0.6
        prop%sigma=0.060
        prop%mw=0.018
        prop%rgas=461
        prop%pvap_coef = [8.1,1770.,240.]
      case(WAT_VAP)
        prop%rho=0.6
        prop%rho0=prop%rho
        prop%mu=1e-5
        prop%cp=2000
        prop%tc=0.03
        prop%sigma=0.060
        prop%mw=0.018
        prop%rgas=461
        prop%pvap_coef = [8.1,1770.,240.]
    end select
    return

    allocate(prop%diff(ne))

    select case(icomponent)
      case(WAT_LIQ)
        prop%rgas=1.0
        prop%mw=1.0
        prop%pvap_coef = (/1.0,0.0,0.0/)
        prop%rho_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%mu_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%cp_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%tc_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%diff_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%sigma_coef = (/1.0,0.0,0.0,0.0,0.0/)
       case(WAT_VAP)
        prop%rgas=1.0
        prop%mw=1.0
        prop%pvap_coef = (/1.0,0.0,0.0/)
        prop%rho_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%mu_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%cp_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%tc_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%diff_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%sigma_coef = (/1.0,0.0,0.0,0.0,0.0/)
       case(AIR)
        prop%rgas=1.0
        prop%mw=1.0
        prop%pvap_coef = (/1.0,0.0,0.0/)
        prop%rho_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%mu_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%cp_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%tc_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%diff_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%sigma_coef = (/1.0,0.0,0.0,0.0,0.0/)
      case default
        prop%rgas=1.0
        prop%mw=1.0
        prop%pvap_coef = (/1.0,0.0,0.0/)
        prop%rho_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%mu_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%cp_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%tc_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%diff_coef = (/1.0,0.0,0.0,0.0,0.0/)
        prop%sigma_coef = (/1.0,0.0,0.0,0.0,0.0/)
    endselect

    end subroutine init_properties

    subroutine calc_pvap(pvap,coef,T,ne)
      real, dimension(ne) :: T,pvap
      real :: coef(3)
      integer :: ne,e
      ! Antoine equation

      do e=1,ne
        pvap(e)= 10.**(coef(1)-coef(2)/(coef(3)+T(e)) )
      end do

    end subroutine

    function poly4(coef,x) result(y)
      real :: x,y
      real :: coef(5),ref

      y=coef(1)+coef(2)*x+coef(3)*x**2+coef(4)*x**3+coef(5)*x**4

    end function

    subroutine calc_poly4(y,coef,x,ne)
      real, dimension(ne) :: x,y
      real :: coef(5)
      integer :: ne,e

      do e=1,ne
        y(e)=poly4(coef,x(e))
      end do
    end subroutine

    subroutine integrate_poly4(y,coef,x,ne,ref)
      real, dimension(ne) :: x,y
      real :: coef(5),ref
      integer :: ne,e,y0
      ! y=int(a+bx+cx^2+dx^3+ex^4)
      y0=coef(1)*ref+coef(2)/2.*ref**2+coef(3)/3.*ref**3+coef(4)/4.*ref**4+coef(5)/5.*ref**5
      do e=1,ne
        y(e)=coef(1)*x(e)+coef(2)/2.*x(e)**2+coef(3)/3.*x(e)**3+coef(4)/4.*x(e)**4+coef(5)/5.*x(e)**5 - y0
      end do

    end subroutine

    subroutine calc_rho_idea(rho,Rgas,P,T,ne)
      real, dimension(ne) :: T,P,rho
      real :: Rgas
      integer :: ne,e

      do e=1,ne
        rho(e)=P(e)/Rgas*T(e)
      end do
    end subroutine

    subroutine calc_temperature(T,energy,cp,ne)
      integer :: ne,e
      real :: energy(*),T(*),cp(*)

      do e=1,ne
        T(e)=energy(e)/cp(e)
      end do

    end subroutine

end module mod_properties
