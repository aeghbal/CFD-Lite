module mod_vfr
  use mod_eqn_setup
  implicit none

  type, extends(equation_t) :: vfr_t
     real, pointer, dimension(:) :: uip,mip
     real, allocatable, dimension(:) :: mfr
     real :: vfr_bg=1e-10
   contains
    procedure :: calc_coef_vfr
  end type

contains
! init vfr field
  function construct_vfr(geom,mip,uip,prop,phase) result(eqn)
    implicit none
    type(vfr_t), pointer :: eqn
    type(properties_t) :: prop
    real, target, dimension(:) :: mip,uip
    type(geometry_t) :: geom
    integer :: ne,nbf,length,phase
    type(bc_t), pointer :: bcp
    character(len=8) :: string

    write(string,'(I8)') phase
    allocate(eqn)
    length = geom%ne+geom%nbf
    eqn%name='vfr_P'//trim(adjustl(string))
    allocate(eqn%phi(length))
    allocate(eqn%phi0(length))
    allocate(eqn%grad(3*length))
    allocate(eqn%mfr(length))
    eqn%uip=>uip
    eqn%mip=>mip
! initialize
    if(phase==GAS) then
      eqn%phi=1.-eqn%vfr_bg
      eqn%mfr=1.
      eqn%phi0=eqn%phi
    else
      eqn%phi=eqn%vfr_bg
      eqn%mfr=0.
      eqn%phi0=eqn%phi
    endif
    eqn%grad=0.

! set BC for vfr
    allocate(eqn%bcs(geom%mg%nintf))
! set BC for vfr
    bcp=>eqn%make_bc(geom,'top');bcp%coef=>zero_flux
    bcp=>eqn%make_bc(geom,'west');bcp%coef=>zero_flux
    bcp=>eqn%make_bc(geom,'east');bcp%coef=>zero_flux
    bcp=>eqn%make_bc(geom,'south');bcp%coef=>zero_flux
    bcp=>eqn%make_bc(geom,'north');bcp%coef=>zero_flux
    bcp=>eqn%make_bc(geom,'bottom');bcp%coef=>zero_flux

  end function
! destroy vfr related coef
  subroutine destroy_vfr(eqn)
    implicit none
    type(vfr_t),pointer :: eqn

    deallocate(eqn%bcs,eqn%grad,eqn%phi0,eqn%phi,eqn%mfr)
    deallocate(eqn)

  end subroutine
! make coef for vfr
  subroutine calc_coef_vfr(eqn,ap,anb,b,geom,prop,dt,src_nucl,src_mass,src_sgn)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(vfr_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn,lf
    real :: dt,src_sgn,phi_ip
    real :: d,f,fnb,sumf,vol,ap(*),anb(*),b(*),sumdefc
    real, dimension(*) :: src_nucl,src_mass
    real :: area,ap0,wt,ds
    real, dimension(3) :: dr,norm,rip,rp,rpnb,gip
    integer :: i,ibc
    real :: beta,dphi,dphic,rd

    do e=1,geom%ne

      ap(e)=0.
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
        phi_ip=(1.-wt)*eqn%phi(e)+wt*eqn%phi(enb)

        ! advection term
        ! inward flux
        f=-fg_sgn*eqn%mip(fg)
        ! upwind bias
        fnb=max(f,0.)

        gip=(1.-wt)*eqn%grad(3*e-2:3*e)+wt*eqn%grad(3*enb-2:3*enb)
        dphi=rsgn(f)*dot_product(gip,-dr)
        dphic=rsgn(f)*(eqn%phi(e)-eqn%phi(enb))
        if((getexpnt(phi_ip)-getexpnt(dphic))<12 .and. (getexpnt(dphic)-getexpnt(tiny(f)))>1.) then
          rd=2*dphi/dphic-1.
          beta=(rd+abs(rd))/(1.+abs(rd))
          sumdefc=sumdefc+0.5*f*beta*dphic
        endif

        anb(idx)=fnb
        ap(e)=ap(e)+fnb
      enddo

      sumf=0.
!      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
!        fg=geom%ef2nb(idx,2)
!        fg_sgn = sgn(fg)
!        fg=abs(fg)
!        sumf=sumf-fg_sgn*eqn%mip(fg)
!      enddo

      ap0=prop%rho(e)*geom%vol(e)/dt
      !ap0=geom%vol(e)/dt
      ap(e)=ap(e)+ap0
      b(e)=ap0*eqn%phi0(e)+sumf*eqn%phi(e)+sumdefc+src_sgn*(src_nucl(e)+src_mass(e))

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
            d=0!prop%tc(e)*area/ds/prop%cp(e)

            !b(e)=b(e)
          case('zero_flux')
            f=0.
            d=0.
            ! TODO
            !d=prop%tc(e)*area/ds/prop%cp(e)
        end select
        ap(e)=ap(e)+d+f
        anb(idx)=anb(idx)+d+f
      end do
    end do
  end subroutine

! boundary conditions for vfr

    subroutine zero_flux(bc,geom,eqn,prop)
      use mod_util

      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='zero_flux'
      select type(eqn)
      type is(vfr_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%phi(e)=eqn%phi(enb)
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
      type is(vfr_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%phi(e)=0.
        end do
      end select

    end subroutine

end module
