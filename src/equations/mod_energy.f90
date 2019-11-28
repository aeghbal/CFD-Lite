module mod_energy
  use mod_eqn_setup
  implicit none

  type, extends(equation_t) :: energy_t
     real, allocatable, dimension(:) :: t,gt
     real, pointer, dimension(:) :: mip
   contains
    procedure :: calc_coef_energy
  end type

contains
! init energy field
  function construct_energy(geom,mip,prop) result(eqn)
    implicit none
    type(energy_t), pointer :: eqn
    type(properties_t) :: prop
    real, target, dimension(:) :: mip
    type(geometry_t) :: geom
    integer :: ne,nbf,length
    type(bc_t), pointer :: bcp

    allocate(eqn)
    length = geom%ne+geom%nbf
    eqn%name='energy'
    allocate(eqn%phi(length))
    allocate(eqn%phi0(length))
    allocate(eqn%grad(3*length))
    allocate(eqn%t(length))
    allocate(eqn%gt(3*length))
    eqn%mip=>mip
! initialize
    eqn%t=273.
    eqn%phi=eqn%t*prop%cp
    eqn%phi0=eqn%phi
    eqn%grad=0.
    eqn%gt=0.
! set BC for energy
    allocate(eqn%bcs(geom%mg%nintf_c2b))
! set BC for energy
    bcp=>eqn%make_bc(geom,'top');bcp%coef=>lid
    bcp=>eqn%make_bc(geom,'west');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'east');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'south');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'north');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'bottom');bcp%coef=>dirichlet0

  end function
! destroy energy related coef
  subroutine destroy_energy(eqn)
    implicit none
    type(energy_t),pointer :: eqn

    deallocate(eqn%bcs,eqn%grad,eqn%phi0,eqn%phi,eqn%t,eqn%gt)
    deallocate(eqn)

  end subroutine

    subroutine solve_energy(eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(energy_t) :: eqn
    type(properties_t) :: prop
    type(geometry_t) :: geom
    real :: dt
    real, dimension(*) :: ap,anb,b,phic
    integer :: nit

    call calc_grad(eqn%t,eqn%gt,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_energy(eqn,ap,anb,b,geom,prop,dt)

    call solve_gs('e',eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    call calc_temperature(eqn%t,eqn%phi,prop%cp,geom%ne)

  end subroutine
! make coef for energy
  subroutine calc_coef_energy(eqn,ap,anb,b,geom,prop,dt)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(energy_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn,lf
    real :: dt
    real :: d,f,fnb,sumf,vol,ap(*),anb(*),b(*),sumdefc
    real :: area,ap0,wt,tci,cpi,ds
    real, dimension(3) :: dr,norm,rip,rp,rpnb,ghi,gti
    integer :: i,ibc
    ! dPHi/dt+Div.(u.Phi)=Div.(dcoef.Grad(Phi))

    do e=1,geom%ne

      ap(e)=0._8
      sumf=0.
      rp=[geom%xc(e),geom%yc(e),geom%zc(e)]
      sumdefc=0.
      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
        anb(idx)=0.
        call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)
        if(lfnb==0) cycle
        fg=geom%ef2nb(idx,2)
        fg_sgn = sgn(fg)
        fg=abs(fg)
        i=3*fg-2
        area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
        norm=fg_sgn*geom%aip(i:i+2)/area
        rip=geom%rip(i:i+2)
        rpnb=[geom%xc(enb),geom%yc(enb),geom%zc(enb)]
        dr = rpnb-rp
        call vec_weight(wt,rip,rp,rpnb)

        ! advection term
        ! inward flux
        f=-fg_sgn*eqn%mip(fg)
        ! upwind bias
        fnb=max(f,0.)
        sumf=sumf+f
        ! diffusion term
        tci=(1.-wt)*prop%tc(e)+wt*prop%tc(enb)
        cpi=(1.-wt)*prop%cp(e)+wt*prop%cp(enb)
        d=tci/cpi/dot_product(dr,norm)*area
        ghi=(1.-wt)*eqn%grad(3*e-2:3*e)+wt*eqn%grad(3*enb-2:3*enb)
        gti=(1.-wt)*eqn%gt(3*e-2:3*e)+wt*eqn%gt(3*enb-2:3*enb)
        sumdefc=sumdefc+tci*area*(dot_product(gti,norm)-dot_product(ghi,norm)/cpi)

        anb(idx)=d+fnb
        ap(e)=ap(e)+d+fnb
      enddo

      ap0=prop%rho(e)*geom%vol(e)/dt
      ap(e)=ap(e)+ap0
      b(e)=ap0*eqn%phi0(e)+sumf*eqn%phi(e)+sumdefc

    end do

   ! set boundary conditions
    do ibc=1,size(eqn%bcs)
      do enb=eqn%bcs(ibc)%esec(1),eqn%bcs(ibc)%esec(2)
        call get_idx(abs(geom%mg%fine_lvl%bs(enb)),0,e,lf)
        idx =geom%ef2nb_idx(e)+lf-1
        fg = geom%ef2nb(idx,2)! always outward on boundary
        i=3*fg-2
        area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
        norm=geom%aip(i:i+2)/area
        rip=geom%rip(i:i+2)
        dr=[geom%xc(enb)-geom%xc(e),geom%yc(enb)-geom%yc(e),geom%zc(enb)-geom%zc(e)]
        ds=dot_product(dr,norm)
        select case(trim(eqn%bcs(ibc)%bc_type))
          case('dirichlet')
            f=0.
            d=prop%tc(e)*area/ds/prop%cp(e)

            !b(e)=b(e)
          case('zero_flux')
            f=0.
            ! TODO
            !d=prop%tc(e)*area/ds/prop%cp(e)
        end select
        ap(e)=ap(e)+d+f
        anb(idx)=anb(idx)+d+f
      end do
    end do
  end subroutine

! boundary conditions for energy

    subroutine lid(bc,geom,eqn,prop)
      use mod_util

      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='dirichlet'
      select type(eqn)
      type is(energy_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%t(e)=373.
          eqn%phi(e)=prop%cp(enb)*eqn%t(e)
        end do
      end select

    end subroutine

    subroutine dirichlet0(bc,geom,eqn,prop)
      use mod_util
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='dirichlet'
      select type(eqn)
      type is(energy_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%t(e)=273.
          eqn%phi(e)=prop%cp(enb)*eqn%t(e)
        end do
      end select

    end subroutine

end module
