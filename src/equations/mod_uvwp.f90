module mod_uvwp
  use mod_eqn_setup
  use mod_util
  implicit none

  type, extends(equation_t) :: uvwp_t
    ! put any extra equation specific variables here
    real, allocatable, dimension(:) :: u,v,w,p
    real, allocatable, dimension(:) :: gu,gv,gw,gp,gpc
    real, allocatable, dimension(:) :: mip,mip0
    real, allocatable, dimension(:) :: u0,v0,w0
    real, allocatable, dimension(:) :: bu,bv,bw,d,dc

   contains
    procedure :: calc_coef_uvw
  end type

contains
! init uvwp field
  function construct_uvwp(geom,prop,dt) result(eqn)
    use mod_properties
    implicit none
    type(uvwp_t), pointer :: eqn
    type(properties_t) :: prop
    type(geometry_t) :: geom
    integer :: ne,nbf,length
    type(bc_t), pointer :: bcp
    real :: dt

    allocate(eqn)
    length = geom%ne+geom%nbf
    eqn%name='uvwp'
    allocate(eqn%u(length))
    allocate(eqn%v(length))
    allocate(eqn%w(length))
    allocate(eqn%u0(length))
    allocate(eqn%v0(length))
    allocate(eqn%w0(length))
    allocate(eqn%p(length))
    !allocate(eqn%guvw(9*length))
    allocate(eqn%gu(3*length))
    allocate(eqn%gv(3*length))
    allocate(eqn%gw(3*length))

    allocate(eqn%gp(3*length))
    allocate(eqn%gpc(3*length))
    length=geom%nf
    allocate(eqn%mip(length))
    allocate(eqn%mip0(length))
    length=geom%ne
    allocate(eqn%bu(length))
    allocate(eqn%bv(length))
    allocate(eqn%bw(length))
    allocate(eqn%d(length))
    allocate(eqn%dc(length))

! initialize
    eqn%u=0.
    eqn%v=0.
    eqn%w=0.
    eqn%p=0.
    eqn%u0=0.
    eqn%v0=0.
    eqn%w0=0.
    eqn%gp=0.
    eqn%gu=0.
    eqn%gv=0.
    eqn%gw=0.
    eqn%mip=0.

    allocate(eqn%bcs(geom%mg%nintf_c2b))
! set BC for uvwp
    bcp=>eqn%make_bc(geom,'top');bcp%coef=>lid
    bcp=>eqn%make_bc(geom,'west');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'east');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'south');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'north');bcp%coef=>dirichlet0
    bcp=>eqn%make_bc(geom,'bottom');bcp%coef=>dirichlet0

    !init mip
    call calc_mip(eqn,prop,geom,.false.,dt)
    eqn%mip0=eqn%mip

  end function
! destroy uvwp related coef
  subroutine destroy_uvwp(eqn)
    implicit none
    type(uvwp_t),pointer :: eqn

    deallocate(eqn%bcs,eqn%u,eqn%v,eqn%w,eqn%p,eqn%u0,eqn%v0,eqn%w0,eqn%gp,eqn%gpc,eqn%mip,eqn%bu,eqn%bv,eqn%bw,eqn%d,eqn%dc)
    deallocate(eqn)

  end subroutine

  subroutine solve_uvwp(eqn,prop,geom,dt,nit,ap,anb,b,phic,subdomain,intf,nsubd)
    use mod_properties
    use mod_solver
    implicit none

    type(uvwp_t) :: eqn
    type(properties_t) :: prop
    type(geometry_t) :: geom
    real :: dt,pref
    real, dimension(*) :: ap,anb,b,phic
    integer :: nit,nsubd
    type(subdomain_t) :: subdomain(nsubd)
    type(intf_t) :: intf(nsubd,nsubd)

!    call calc_grad(eqn%p,eqn%gp,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_uvw(eqn,ap,anb,geom,prop,dt)

    !call set_a_0(eqn%u,geom%ne);call set_a_0(eqn%v,geom%ne);call set_a_0(eqn%w,geom%ne)
    call solve_gs('u',eqn%u,ap,anb,eqn%bu,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
    call solve_gs('v',eqn%v,ap,anb,eqn%bv,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
    call solve_gs('w',eqn%w,ap,anb,eqn%bw,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    call calc_grad(eqn%u,eqn%gu,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call calc_grad(eqn%v,eqn%gv,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call calc_grad(eqn%w,eqn%gw,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
   ! Rhie Chow
    call calc_mip(eqn,prop,geom,.true.,dt)

    call calc_coef_p(eqn,ap,anb,b,geom,prop)

    call set_a_0(phic,geom%ne+geom%nbf)
    call solve('pc',subdomain,intf,geom,ap,anb,b,phic,nsubd,nit)

    pref=phic(1)! if closed system
    call adjust_pc(phic,pref,geom)
    call calc_grad(phic,eqn%gpc,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call update_uvwp(eqn,prop,geom,phic,geom%ne,geom%nf,geom%nbf)

  end subroutine

  subroutine adjust_pc(pc,pref,geom)
    use mod_solver
    type(geometry_t) :: geom
    integer :: e,enb,lfnb,idx,lf,enb1
    real :: pc(geom%ne+geom%nbf),pref
    real :: ga11,ga12,ga13,ga22,ga23,ga33,A(3,3),B(3,1),C(3,1),A_inv(3,3),dphi,dr(3),wt,grad(3)

    do e=1,geom%ne
      pc(e)=pc(e)-pref
    end do
    do enb=geom%ne+1,geom%ne+geom%nbf
      call get_idx(abs(geom%mg%fine_lvl%bs(enb)),0,e,lf)

      if(.false.) then! first order extrapolation
        dr = [geom%xc(enb)-geom%xc(e),geom%yc(enb)-geom%yc(e),geom%zc(enb)-geom%zc(e)]
        call internal_extrapolation_grad(pc,e,grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
        pc(enb)=pc(e)-dot_product(grad,dr)
      else! zeroth order extrapolation
        pc(enb)=pc(e)
      endif
    end do

  end subroutine

! make coef for uvwp
  subroutine calc_coef_uvw(eqn,ap,anb,geom,prop,dt)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(uvwp_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn,lf
    real :: dt,ap(geom%ne),anb(geom%nf*2-geom%nbf)
    real :: d,f,fnb,sumf,vol,wt,muip,ds,ds_p,area,ap0
    real, dimension(3) :: dr,norm,rip,vrel,vbnc,gip,rp,rpnb,rp_p,rpnb_p,sumss,sumdefc,dr_p,drip
    integer :: i,ibc,m

    do e=1,geom%ne

      rp=[geom%xc(e),geom%yc(e),geom%zc(e)]
      ap(e)=0.
      sumf=0.;sumss=0.;sumdefc=0.
      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
        anb(idx)=0.
        call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)
        f=0.;fnb=0.;d=0.
        if(lfnb>0) then
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

          drip=rip-rp
          rp_p=rip-dot_product(drip,norm)*norm
          drip=rip-rpnb
          rpnb_p=rip-dot_product(drip,norm)*norm
          dr_p=rpnb_p-rp_p
          ds_p=sqrt(dot_product(dr_p,dr_p))
          ! advection term
          ! inward flux
          f=-fg_sgn*eqn%mip(fg)
          ! upwind bias
          fnb=max(f,0.)
          sumf=sumf+f
          ! diffusion term
          muip=(1.-wt)*prop%mu(e)+wt*prop%mu(enb)
          d=muip*area/ds

          ! secondary stress term
          do m=1,3
            gip(1)=(1.-wt)*eqn%gu(3*e-3+m)+wt*eqn%gu(3*enb-3+m)
            gip(2)=(1.-wt)*eqn%gv(3*e-3+m)+wt*eqn%gv(3*enb-3+m)
            gip(3)=(1.-wt)*eqn%gw(3*e-3+m)+wt*eqn%gw(3*enb-3+m)
            sumss(m)=sumss(m)+muip*area*dot_product(gip,dr)/ds
          enddo
          ! Deferred correction of real diffusion
          gip=(1.-wt)*eqn%gu(3*e-2:3*e)+wt*eqn%gu(3*enb-2:3*enb)
          sumdefc(1)=sumdefc(1)+muip*area*(dot_product(gip,dr_p)/ds_p-dot_product(gip,dr)/ds)
          gip=(1.-wt)*eqn%gv(3*e-2:3*e)+wt*eqn%gv(3*enb-2:3*enb)
          sumdefc(2)=sumdefc(2)+muip*area*(dot_product(gip,dr_p)/ds_p-dot_product(gip,dr)/ds)
          gip=(1.-wt)*eqn%gw(3*e-2:3*e)+wt*eqn%gw(3*enb-2:3*enb)
          sumdefc(3)=sumdefc(3)+muip*area*(dot_product(gip,dr_p)/ds_p-dot_product(gip,dr)/ds)

        endif

        anb(idx)=d+fnb
        ap(e)=ap(e)+d+fnb
      enddo

      ap0=prop%rho(e)*geom%vol(e)/dt
      ap(e)=ap(e)+ap0

      eqn%bu(e)=ap0*eqn%u0(e)+sumf*eqn%u(e)-geom%vol(e)*eqn%gp(3*e-2)+sumss(1)+sumdefc(1)
      eqn%bv(e)=ap0*eqn%v0(e)+sumf*eqn%v(e)-geom%vol(e)*eqn%gp(3*e-1)+sumss(2)+sumdefc(2)
      eqn%bw(e)=ap0*eqn%w0(e)+sumf*eqn%w(e)-geom%vol(e)*eqn%gp(3*e)  +sumss(3)+sumdefc(3)

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
        ds=sqrt(dot_product(dr,dr))
        select case(trim(eqn%bcs(ibc)%bc_type))
          case('dirichlet')
            f=0.
            d=prop%mu(e)*area/ds
            vbnc=[eqn%u(enb),eqn%v(enb),eqn%w(enb)]
            vrel=[eqn%u(e),eqn%v(e),eqn%w(e)]
            vrel=vrel-dot_product(vrel,norm)*norm! transverse velocity
            vrel=vbnc-vrel

            eqn%bu(e)=eqn%bu(e)+d*vrel(1)-d*eqn%u(e)
            eqn%bv(e)=eqn%bv(e)+d*vrel(2)-d*eqn%v(e)
            eqn%bw(e)=eqn%bw(e)+d*vrel(3)-d*eqn%w(e)
          case('zero_flux')
            f=0.
            d=prop%mu(e)*area/ds
        end select
        ap(e)=ap(e)+d+f
        anb(idx)=anb(idx)+d+f
      end do
    end do

    ! calculate Rhie-Chow and SIMPLEC coef
    do e=1,geom%ne
      eqn%d(e)=geom%vol(e)/ap(e)! Rhie-Chow

      eqn%dc(e)=ap(e)
      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
        eqn%dc(e)=eqn%dc(e)-anb(idx)
      enddo
      eqn%dc(e)=geom%vol(e)/eqn%dc(e)! SIMPLEC
    enddo

  end subroutine

! make coef for uvwp
  subroutine calc_coef_p(eqn,ap,anb,b,geom,prop)
    use mod_util
    use mod_properties
    implicit none
    type(properties_t) :: prop
    type(geometry_t) :: geom
    class(uvwp_t) :: eqn
    integer :: e,enb,lfnb,idx,fg,fg_sgn,lf,idxnb
    real, dimension(geom%ne) :: ap,b
    real, dimension(2*geom%nf-geom%nbf) :: anb
    real :: dt,wt
    real :: d,f,fnb,sumf,vol,area,rhoip
    real, dimension(3) :: rp,rpnb,dr,norm,rip
    integer :: i,ibc


    do e=1,geom%ne

      ap(e)=0.
      sumf=0.
      rp=[geom%xc(e),geom%yc(e),geom%zc(e)]

      do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
        anb(idx)=0.
        call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)

        d=0.
        if(lfnb>0) then
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

          ! inward flux
          f=-fg_sgn*eqn%mip(fg)
          sumf=sumf+f

          ! diffusion term
          rhoip=(1.-wt)*prop%rho(e)+wt*prop%rho(enb)
          d=((1.-wt)*eqn%dc(e)+wt*eqn%dc(enb))/dot_product(dr,norm)*rhoip*area
        endif

        anb(idx)=d
        ap(e)=ap(e)+d
      enddo
      b(e)=sumf

    end do
    ! set coefficients for boundary condition
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
        select case(trim(eqn%bcs(ibc)%bc_type))
          case('dirichlet')
            d=0.

          case('zero_flux')
            d=0.

        end select
        ap(e)=ap(e)+d
        anb(idx)=anb(idx)+d
        b(e)=b(e)-eqn%mip(fg)
      end do
    end do

  end subroutine

  subroutine update_uvwp(eqn,prop,geom,pc,ne,nf,nbf)
    use mod_properties
    use mod_util
    implicit none
    type(properties_t) :: prop
    type(uvwp_t) :: eqn
    type(geometry_t) :: geom
    integer :: ne,nf,nbf
    real, dimension(ne+nbf) :: pc
    real :: dmip,rhoip,area,dip,wt
    real, dimension(3) :: norm,rp,rpnb,dr,rip
    integer :: e,fg,idx,enb,lfnb,lf,i,ibc,idx1

    do e=1,ne
      eqn%p(e)= eqn%p(e)+pc(e)
      eqn%gp(3*e-2:3*e)=eqn%gp(3*e-2:3*e)+eqn%gpc(3*e-2:3*e)
      ! correct velocity field
      if(.false.) then
      eqn%u(e)=eqn%u(e)-eqn%dc(e)*eqn%gpc(3*e-2)
      eqn%v(e)=eqn%v(e)-eqn%dc(e)*eqn%gpc(3*e-1)
      eqn%w(e)=eqn%w(e)-eqn%dc(e)*eqn%gpc(3*e)
      endif
    end do

    do fg=1,nf
      call get_idx(geom%mg%fine_lvl%s2g(fg),0,e,lf)
      idx=geom%ef2nb_idx(e)+lf-1
      call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)
      if(lfnb==0) cycle

      i=3*fg-2
      area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
      norm=geom%aip(i:i+2)/area
      rp=[geom%xc(e),geom%yc(e),geom%zc(e)]
      rpnb=[geom%xc(enb),geom%yc(enb),geom%zc(enb)]
      dr=rpnb-rp
      rip=geom%rip(i:i+2)
      call vec_weight(wt,rip,rp,rpnb)

      dmip=0.
      dip=(1.-wt)*eqn%dc(e)+wt*eqn%dc(enb)
      rhoip=(prop%rho(e)+prop%rho(enb))/2.
      dmip=rhoip*area*dip*(pc(enb)-pc(e))/dot_product(dr,norm)
      eqn%mip(fg)=eqn%mip(fg)-dmip

    end do

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
        select case(trim(eqn%bcs(ibc)%bc_type))
          case('dirichlet')

          case('zero_flux')

        end select
      end do
    end do

  end subroutine

  subroutine calc_mip(eqn,prop,geom,lRhieChow,dt)
    use mod_util
    use mod_properties
    implicit none
    type(uvwp_t) :: eqn
    type(geometry_t) :: geom
    type(properties_t) :: prop
    real :: wt,area,veln,dip,rhoip,dt
    real, dimension(3) :: rp,rpnb,rip,norm,vec1,vec2,dr,velip,gpip,velip0
    integer :: fg,e,lf,idx,enb,lfnb,i
    logical :: lRhieChow

    associate(mip=>eqn%mip,u=>eqn%u,v=>eqn%v,w=>eqn%w)

    do fg=1,geom%nf
      call get_idx(geom%mg%fine_lvl%s2g(fg),0,e,lf)
      idx=geom%ef2nb_idx(e)+lf-1
      call get_idx(geom%ef2nb(idx,1),0,enb,lfnb)
      if(lfnb==0) cycle

      rp=[geom%xc(e),geom%yc(e),geom%zc(e)]
      rpnb=[geom%xc(enb),geom%yc(enb),geom%zc(enb)]
      i=3*fg-2
      area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
      norm=geom%aip(i:i+2)/area
      rip=geom%rip(i:i+2)
      call vec_weight(wt,rip,rp,rpnb)
      !call vec_weight_proj(wt,rip,rp,rpnb,norm)

      vec1=[u(e),v(e),w(e)]
      vec2=[u(enb),v(enb),w(enb)]
      velip= (1.-wt)*vec1+wt*vec2
      rhoip=prop%rho(e)*(1.-wt)+prop%rho(enb)*wt

      mip(fg)=dot_product(velip,norm)*rhoip*area

      ! Rhie-Chow
      if(lRhieChow) then ! internal faces
        dr=rpnb-rp
        gpip =(1.-wt)*eqn%gp(3*e-2:3*e)+wt*eqn%gp(3*enb-2:3*enb)
        dip=(1.-wt)*eqn%d(e)+wt*eqn%d(enb)
        vec1=[eqn%u0(e),eqn%v0(e),eqn%w0(e)]
        vec2=[eqn%u0(enb),eqn%v0(enb),eqn%w0(enb)]
        velip0= (1.-wt)*vec1+wt*vec2

        mip(fg) = mip(fg)-rhoip*area*dip/dot_product(dr,norm)*(eqn%p(enb) - eqn%p(e) - dot_product(gpip,dr))  &
                         -rhoip/dt*dip*(eqn%mip0(fg)-dot_product(velip0,norm)*rhoip*area)
      endif

    end do
    end associate

  end subroutine

! boundary conditions for energy
    subroutine dirichlet0(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='dirichlet'
      select type(eqn)
      type is(uvwp_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%u(e)=0.
          eqn%v(e)=0.
          eqn%w(e)=0.
          eqn%p(e)=eqn%p(enb)
          idx =geom%ef2nb_idx(enb)+lfnb-1
          fg = geom%ef2nb(idx,2)
          eqn%mip(fg)=0.
        end do
      end select

    end subroutine

    subroutine symmetry(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,idx,enb,lf,lfnb,fg,i
      real :: area, norm(3),vel(3),veln(3),velt(3)

      bc%bc_type='zero_flux'
      select type(eqn)
      type is(uvwp_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          idx =geom%ef2nb_idx(enb)+lfnb-1
          fg = geom%ef2nb(idx,2)
          i=3*fg-2
          area=sqrt(geom%aip(i)**2+geom%aip(i+1)**2+geom%aip(i+2)**2)
          norm=geom%aip(i:i+2)/area
          vel=[eqn%u(enb),eqn%v(enb),eqn%w(enb)]
          veln=dot_product(vel,norm)*norm
          velt=vel-veln
          eqn%u(e)=velt(1)-2*veln(1)
          eqn%v(e)=velt(2)-2*veln(2)
          eqn%w(e)=velt(3)-2*veln(3)
          eqn%p(e)=eqn%p(enb)

          eqn%mip(fg)=0.
        end do
      end select
    end subroutine

    subroutine lid(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e,enb,lfnb,idx,fg

      bc%bc_type='dirichlet'
      select type(eqn)
      type is(uvwp_t)
        do e=bc%esec(1),bc%esec(2)
          call get_idx(abs(geom%mg%fine_lvl%bs(e)),0,enb,lfnb)
          eqn%u(e)=1.
          eqn%v(e)=0.
          eqn%w(e)=0.
          eqn%p(e)=eqn%p(enb)
          idx =geom%ef2nb_idx(enb)+lfnb-1
          fg = geom%ef2nb(idx,2)
          eqn%mip(fg)=0.
        end do
      end select

    end subroutine

    subroutine inlet(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e

      bc%bc_type='inlet'
      select type(eqn)
      type is(uvwp_t)
        do e=bc%esec(1),bc%esec(2)
          eqn%u(e)=0.
          eqn%v(e)=0.
          eqn%w(e)=-1.
        end do
      end select

    end subroutine

    subroutine outlet(bc,geom,eqn,prop)
      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop
      integer :: e

      bc%bc_type='outlet_p'
      select type(eqn)
      type is(uvwp_t)
        do e=bc%esec(1),bc%esec(2)
          eqn%p(e)=0.
        end do
      end select

    end subroutine


end module
