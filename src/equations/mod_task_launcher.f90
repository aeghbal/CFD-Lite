module mod_task_launcher
  use mod_physics
  use mod_cell
  implicit none
contains
   subroutine solve_uvw(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_cell
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,n,m

  do n=1,NPHASES
    phys%iphase=n
    associate(eqn=>phys%phase(n)%uvwp,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt,alpha=>phys%phase(n)%vfr%phi)

    ! update boundaries
    do m=1,phys%nintf
      call eqn%bcs(m)%coef(geom,eqn,prop)
    enddo

    if(phys%iphase==GAS) then
  !    call calc_grad(eqn%p,eqn%gp,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

      call calc_coef_uvw(eqn,ap,anb,geom,prop,dt, &
                        phys%interaction(LIQUID,GAS)%scheme(MASS_EXCH)%src,phys%interaction(LIQUID,GAS)%scheme(MOMENTUM_EXCH)%sdc,phys%interaction(LIQUID,GAS)%scheme(MASS_EXCH)%sgn,alpha)

      !call set_a_0(eqn%u,geom%ne);call set_a_0(eqn%v,geom%ne);call set_a_0(eqn%w,geom%ne)
      call solve_gs('u',eqn%u,ap,anb,eqn%bu,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
      call solve_gs('v',eqn%v,ap,anb,eqn%bv,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
      call solve_gs('w',eqn%w,ap,anb,eqn%bw,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

      call calc_grad(eqn%u,eqn%gu,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call calc_grad(eqn%v,eqn%gv,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call calc_grad(eqn%w,eqn%gw,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
   ! Rhie Chow
      call calc_mip(eqn,prop,geom,.true.,dt,phys%phase(n)%vfr%phi)

    elseif(phys%iphase==LIQUID) then

      call calc_coef_uvw(eqn,ap,anb,geom,prop,dt, &
                        phys%interaction(GAS,LIQUID)%scheme(MASS_EXCH)%src,phys%interaction(GAS,LIQUID)%scheme(MOMENTUM_EXCH)%sdc,phys%interaction(GAS,LIQUID)%scheme(MASS_EXCH)%sgn,alpha)

      !call set_a_0(eqn%u,geom%ne);call set_a_0(eqn%v,geom%ne);call set_a_0(eqn%w,geom%ne)
      call solve_gs('u',eqn%u,ap,anb,eqn%bu,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
      call solve_gs('v',eqn%v,ap,anb,eqn%bv,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
      call solve_gs('w',eqn%w,ap,anb,eqn%bw,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

      call calc_grad(eqn%u,eqn%gu,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call calc_grad(eqn%v,eqn%gv,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call calc_grad(eqn%w,eqn%gw,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
   ! Rhie Chow
      call calc_mip(eqn,prop,geom,.true.,dt,phys%phase(n)%vfr%phi)

    else

      bg=AIR

      call calc_uvw_alg (eqn,geom,prop,dt,phys%phase(bg)%vfr%phi,phys%phase(bg)%ndf%phi, &
                        phys%phase(bg)%uvwp%u,phys%phase(bg)%uvwp%v,phys%phase(bg)%uvwp%w,phys%phase(bg)%prop%rho)

      call calc_mip(eqn,prop,geom,.false.,dt,phys%phase(n)%vfr%phi)

    endif

    end associate
  enddo
    ! calculate mixture momentum in preparation for pressure solution
    call calc_mixture_field(phys%phase(MIXTURE)%uvwp,phys%phase,geom)

    ! calculate source terms due to interaction between phases
    do n=1,NPHASES
      do m=n+1,NPHASES
!        call phys%interaction(n,m)%update(NUCLEATION,geom)
!        call phys%interaction(n,m)%update(MASS_EXCH,geom)
        call phys%interaction(n,m)%update(MOMENTUM_EXCH,geom)
!        call phys%interaction(n,m)%update(ENERGY_EXCH,geom)
      end do
    end do

  end subroutine

  subroutine solve_pressure(phys,geom)
    use mod_properties
    use mod_solver
    implicit none

    type(properties_t) :: prop
    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: dt,pref
    integer :: nit,bg

    associate(eqn=>phys%phase(MIXTURE)%uvwp,prop=>phys%phase(MIXTURE)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

      ! update boundaries
!      do m=1,phys%nintf
!        call eqn%bcs(m)%coef(geom,eqn,prop)
!      enddo

      call calc_coef_p(eqn,ap,anb,b,geom,prop)

      call set_a_0(phic,geom%ne+geom%nbf)
      call solve_gs('pc',phic,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

      pref=phic(1)! if closed system
      call adjust_pc(phic,pref,geom)
      call calc_grad(phic,eqn%gpc,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call update_uvwp(eqn,prop,geom,phic,eqn%gpc)
      call update_mip(eqn,prop,geom,phic)

    end associate

    ! update gas mip as well
    phys%phase(GAS)%uvwp%p=phys%phase(MIXTURE)%uvwp%p
    phys%phase(GAS)%uvwp%gp=phys%phase(MIXTURE)%uvwp%gp
    call update_mip(phys%phase(GAS)%uvwp,phys%phase(GAS)%prop,geom,phys%phic)
    phys%phase(LIQUID)%uvwp%p=phys%phase(MIXTURE)%uvwp%p
    phys%phase(LIQUID)%uvwp%gp=phys%phase(MIXTURE)%uvwp%gp
    call update_mip(phys%phase(LIQUID)%uvwp,phys%phase(LIQUID)%prop,geom,phys%phic)
    ! calculate partial pressure
    ! call calc_pvap
  end subroutine

  subroutine solve_energy(phys,geom)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,n,m

  do n=1,NPHASES
    phys%iphase=n
    associate(eqn=>phys%phase(n)%energy,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt,alpha=>phys%phase(n)%vfr%phi)

! update boundaries
    do m=1,phys%nintf
      call eqn%bcs(m)%coef(geom,eqn,prop)
    enddo

    call calc_grad(eqn%t,eqn%gt,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    if(n==GAS) then
      call calc_coef_energy(eqn,ap,anb,b,geom,prop,dt,&
                            phys%interaction(LIQUID,GAS)%scheme(MASS_EXCH)%src,phys%interaction(LIQUID,GAS)%scheme(ENERGY_EXCH)%src,phys%interaction(LIQUID,GAS)%scheme(MASS_EXCH)%sgn,alpha)
    elseif(n==LIQUID) then
      call calc_coef_energy(eqn,ap,anb,b,geom,prop,dt,&
                            phys%interaction(GAS,LIQUID)%scheme(MASS_EXCH)%src,phys%interaction(GAS,LIQUID)%scheme(ENERGY_EXCH)%src,phys%interaction(GAS,LIQUID)%scheme(MASS_EXCH)%sgn,alpha)
    endif
    call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    call calc_temperature(eqn%t,eqn%phi,prop%cp,geom%ne)

    end associate
  enddo
    ! calculate source terms due to interaction between phases
    do n=1,NPHASES
      do m=n+1,NPHASES
!        call phys%interaction(n,m)%update(NUCLEATION,geom)
!        call phys%interaction(n,m)%update(MASS_EXCH,geom)
!        call phys%interaction(n,m)%update(MOMENTUM_EXCH,geom)
        call phys%interaction(n,m)%update(ENERGY_EXCH,geom)
      end do
    end do

  end subroutine

  subroutine solve_vfr(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,n,m

    ! solve phase continuity equations
    do n=1,NPHASES
      phys%iphase=n
      associate(eqn=>phys%phase(n)%vfr,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt,alpha=>phys%phase(n)%vfr%phi)
! update boundaries
      do m=1,phys%nintf
        call eqn%bcs(m)%coef(geom,eqn,prop)
      enddo

      if(n==GAS) then

        call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

        call calc_coef_vfr(eqn,ap,anb,b,geom,prop,dt,phys%interaction(LIQUID,GAS)%scheme(NUCLEATION)%src,phys%interaction(LIQUID,GAS)%scheme(MASS_EXCH)%src,phys%interaction(LIQUID,GAS)%scheme(NUCLEATION)%sgn)

        call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

        call rescale_vfr(eqn%phi,geom%ne)

      elseif(n==LIQUID) then

        call calc_algebraic_vfr(phys%phase,n,geom%ne)

      endif
      end associate
    enddo

    call calc_phase_mfr(phys%phase,geom)

    ! calculate source terms due to interaction between phases
    do n=1,NPHASES
      do m=n+1,NPHASES
!        call phys%interaction(n,m)%update(NUCLEATION,geom)
        call phys%interaction(n,m)%update(MASS_EXCH,geom)
!        call phys%interaction(n,m)%update(MOMENTUM_EXCH,geom)
!        call phys%interaction(n,m)%update(ENERGY_EXCH,geom)
      end do
    end do

  end subroutine

  subroutine solve_mfr(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,p,c,m

    p=phys%iphase
    do c=1,NCOMPONENTS
      if(.not. phys%phase(p)%COMPONENT_SET(c)) cycle
      associate(eqn=>phys%phase(p)%component(c)%mfr,prop=>phys%phase(p)%component(c)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)
! update boundaries
        do m=1,phys%nintf
          call eqn%bcs(m)%coef(geom,eqn,prop)
        enddo

        call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

        call calc_coef_mfr(eqn,ap,anb,b,geom,prop,dt)

        call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

        !call calc_component_vfr(eqn%vfr,eqn%phi,geom,prop,phase%prop)
      end associate
    enddo


  end subroutine


  subroutine solve_scalar(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,n,m

    n=phys%iphase
    associate(eqn=>phys%scalar,prop=>phys%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt,alpha=>phys%phase(n)%vfr%phi)
! update boundaries
    do m=1,phys%nintf
      call eqn%bcs(m)%coef(geom,eqn,prop)
    enddo

    call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_scalar(eqn,ap,anb,b,geom,prop,dt)

    call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    end associate
  end subroutine

  subroutine solve_ndf(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,n,m


    phys%iphase=LIQUID
    n=LIQUID
    associate(eqn=>phys%phase(n)%ndf,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt,alpha=>phys%phase(n)%vfr%phi)

! update boundaries
    do m=1,phys%nintf
      call eqn%bcs(m)%coef(geom,eqn,prop)
    enddo
    !call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_ndf(eqn,ap,anb,b,geom,prop,dt,phys%interaction(GAS,LIQUID)%scheme(NUCLEATION)%src,phys%interaction(GAS,LIQUID)%scheme(NUCLEATION)%sgn)

    call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    end associate


    ! calculate source terms due to interaction between phases
    do n=1,NPHASES
      do m=n+1,NPHASES
        call phys%interaction(n,m)%update(NUCLEATION,geom)
!        call phys%interaction(n,m)%update(MASS_EXCH,geom)
!        call phys%interaction(n,m)%update(MOMENTUM_EXCH,geom)
!        call phys%interaction(n,m)%update(ENERGY_EXCH,geom)
      end do
    end do

  end subroutine

  subroutine solve_properties(phys,geom)!prop,iphase,icomponent,T,P,ne)
    implicit none
    type(phys_t) :: phys
    type(geometry_t) :: geom
    integer :: iphase,icomponent,ne,n


    call update_mixture_properties(phys%phase,geom%ne)
    associate(prop=>phys%phase(GAS)%prop)
      call calc_pvap(prop%pvap,prop%pvap_coef,phys%phase(GAS)%energy%t,ne)
    end associate
    return

    ne=geom%ne
  do n=1,NPHASES
    phys%iphase=n
    associate(prop=>phys%phase(n)%prop)

    ! calculate density
    if(iphase==GAS) then
      call calc_pvap(prop%pvap,prop%pvap_coef,phys%phase(n)%energy%t,ne)
      call calc_rho_idea(prop%rho,prop%rgas,phys%phase(MIXTURE)%uvwp%p,phys%phase(n)%energy%t,ne)
    else
      call calc_poly4(prop%rho,prop%rho_coef,phys%phase(n)%energy%t,ne)
    endif
    call calc_poly4(prop%mu,prop%mu_coef,phys%phase(n)%energy%t,ne)
    call calc_poly4(prop%cp,prop%cp_coef,phys%phase(n)%energy%t,ne)
    call calc_poly4(prop%tc,prop%tc_coef,phys%phase(n)%energy%t,ne)
    call calc_poly4(prop%diff,prop%diff_coef,phys%phase(n)%energy%t,ne)
    call calc_poly4(prop%sigma,prop%sigma_coef,phys%phase(n)%energy%t,ne)

    endassociate
  enddo

  end subroutine
end module
