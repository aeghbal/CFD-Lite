module mod_scalar
  use mod_eqn_setup
  implicit none

  type, extends(equation_t) :: scalar_t
    ! put any extra equation specific variables here
    real :: dcoef=1.
    real :: vel(3)=[0.,0.,-100.]

   contains
    procedure :: calc_coef_scalar
  end type

contains
! init scalar field
  function construct_scalar(geom) result(eqn)
    implicit none
    type(scalar_t), pointer :: eqn
    type(geometry_t) :: geom
    integer :: ne,nbf,length,nintf_c2b
    type(bc_t), pointer :: bcp

    allocate(eqn)
    length = geom%ne+geom%nbf
    eqn%name='scalar'
    allocate(eqn%phi(length))
    allocate(eqn%phi0(length))
    allocate(eqn%grad(3*length))
! initialize
    eqn%phi=0.
    eqn%phi0=0.
    eqn%grad=0.
! set BC for scalar
    nintf_c2b=3
    allocate(eqn%bcs(nintf_c2b))
    bcp=>eqn%make_bc(geom,'fp-connect');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'out_wall');bcp%coef=>dirichlet1
    bcp=>eqn%make_bc(geom,'outlet');bcp%coef=>dirichlet0
    !bcp=>eqn%make_bc(geom,'tri_cc1_inlet');bcp%coef=>dirichlet0
    !bcp=>eqn%make_bc(geom,'quad_cc1_inlet');bcp%coef=>dirichlet0
    !bcp=>eqn%make_bc(geom,'f2sol_walls_CC');bcp%coef=>dirichlet1
    !bcp=>eqn%make_bc(geom,'tri_cc1_outlet');bcp%coef=>dirichlet0
    !bcp=>eqn%make_bc(geom,'quad_cc1_outlet');bcp%coef=>dirichlet0

  end function
! destroy scalar related coef
  subroutine destroy_scalar(eqn)
    implicit none
    type(scalar_t),pointer :: eqn

    deallocate(eqn%bcs,eqn%grad,eqn%phi0,eqn%phi)
    deallocate(eqn)

  end subroutine

  subroutine solve_scalar(eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(scalar_t) :: eqn
    type(properties_t) :: prop
    type(geometry_t) :: geom
    real :: dt
    real, dimension(*) :: ap,anb,b,phic
    integer :: nit
    ! solve scalar
    call eqn%calc_coef_scalar(geom,prop,ap,anb,b,dt)
    ! calc Grad scalar
    call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    ! solve
    call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)


  end subroutine
! make coef for scalar
  subroutine calc_coef_scalar(eqn,geom,prop,ap,anb,b,dt)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(scalar_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn
    real :: dt,wnb
    real :: d,f,fnb,sumf,vol,ap(*),anb(*),b(*)
    real :: dr(3),area,norm(3),rip(3),ap0
    integer :: i
    ! dPHi/dt+Div.(u.Phi)=Div.(dcoef.Grad(Phi))

    do e=1,geom%ne

      ap(e)=0._8
      sumf=0.
      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
        fg=geom%ef2nb(idx,2)
        fg_sgn = sgn(fg)
        fg=abs(fg)
        i=3*fg-2
        area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
        norm=fg_sgn*geom%aip(i:i+2)/area
        rip=geom%rip(i:i+2)

        call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)
        dr=[geom%xc(enb)-geom%xc(e),geom%yc(enb)-geom%yc(e),geom%zc(enb)-geom%zc(e)]
        ! advection term
        ! inward flux
        f=dot_product(eqn%vel,-norm)*area
        ! upwind bias
        wnb=0.
        if(f>0.) wnb=1.
        fnb=wnb*f
        sumf=sumf+f
        ! diffusion term
        d=eqn%dcoef/dot_product(dr,dr)*dot_product(dr,norm)*area

        anb(idx)=d+fnb
        ap(e)=ap(e)+d+fnb
      enddo

      ap0=geom%vol(e)/dt
      ap(e)=ap(e)+ap0
      b(e)=ap0*eqn%phi0(e)+sumf*eqn%phi(e)

    end do


  end subroutine

! boundary conditions for energy
    subroutine dirichlet0(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(properties_t) :: prop
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      integer :: e

      do e=bc%esec(1),bc%esec(2)
        eqn%phi(e)=0.
      end do

    end subroutine

    subroutine dirichlet1(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(properties_t) :: prop
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      integer :: e

      do e=bc%esec(1),bc%esec(2)
        eqn%phi(e)=1.
      end do

    end subroutine

end module

