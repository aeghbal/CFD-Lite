module mod_mfr
  use mod_eqn_setup
  implicit none

  type, extends(equation_t) :: mfr_t
     real, pointer, dimension(:) :: mip
     real, allocatable, dimension(:) :: vfr
   contains
    procedure :: calc_coef_mfr
  end type

contains
! init mfr field
  function construct_mfr(geom,mip,prop) result(eqn)
    implicit none
    type(mfr_t), pointer :: eqn
    type(properties_t) :: prop
    real, target, dimension(:) :: mip
    type(geometry_t) :: geom
    integer :: ne,nbf,length
    type(bc_t), pointer :: bcp

    allocate(eqn)
    length = geom%ne+geom%nbf
    eqn%name='mfr'
    allocate(eqn%phi(length))
    allocate(eqn%phi0(length))
    allocate(eqn%grad(3*length))
    allocate(eqn%vfr(length))
    eqn%mip=>mip
! initialize
    eqn%phi=0.
    eqn%phi0=eqn%phi
    eqn%grad=0.
    eqn%vfr=0.
! set BC for mfr
    allocate(eqn%bcs(geom%mg%nintf))
! set BC for mfr
    bcp=>eqn%make_bc(geom,'top');bcp%coef=>lid
    bcp=>eqn%make_bc(geom,'west');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'east');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'south');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'north');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'bottom');bcp%coef=>dirichlet0

  end function
! destroy mfr related coef
  subroutine destroy_mfr(eqn)
    implicit none
    type(mfr_t),pointer :: eqn

    deallocate(eqn%bcs,eqn%grad,eqn%phi0,eqn%phi)
    deallocate(eqn)

  end subroutine
! make coef for mfr
  subroutine calc_coef_mfr(eqn,ap,anb,b,geom,prop,dt)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(mfr_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn,lf
    real :: dt
    real :: d,f,fnb,sumf,vol,ap(*),anb(*),b(*),sumdefc
    real :: area,ap0,wt,tci,cpi,muip,ds,ds_p
    real, dimension(3) :: dr,norm,rip,rp,rpnb,ghi,drip,rp_p,rpnb_p,dr_p
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
        ds=sqrt(dot_product(dr,dr))
        call vec_weight(wt,rip,rp,rpnb)

        ! advection term
        ! inward flux
        f=-fg_sgn*eqn%mip(fg)
        ! upwind bias
        fnb=max(f,0.)
        sumf=sumf+f
        ! diffusion term
        muip=(1.-wt)*prop%mu(e)+wt*prop%mu(enb)
        d=muip/ds*area

        drip=rip-rp
        rp_p=rip-dot_product(drip,norm)*norm
        drip=rip-rpnb
        rpnb_p=rip-dot_product(drip,norm)*norm
        dr_p=rpnb_p-rp_p
        ds_p=dot_product(dr_p,norm)

        ghi=(1.-wt)*eqn%grad(3*e-2:3*e)+wt*eqn%grad(3*enb-2:3*enb)

        sumdefc=sumdefc+muip*area*(dot_product(ghi,dr_p)/ds_p-dot_product(ghi,dr)/ds)

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
            d=0.!prop%tc(e)*area/ds/prop%cp(e)

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

! boundary conditions for mfr

    subroutine lid(bc,geom,eqn,prop)
      use mod_util

      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='dirichlet'
      select type(eqn)
      type is(mfr_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%phi(e)=0.
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
      type is(mfr_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%phi(e)=0.
        end do
      end select

    end subroutine

end module
