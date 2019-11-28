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
    associate(eqn=>phys%phase(n)%uvwp,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

    if(phys%iphase==GAS) then
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

    elseif(phys%iphase==LIQUID) then

      bg=AIR

      call calc_uvw_alg (eqn,geom,prop,dt,phys%phase(bg)%vfr%phi,phys%phase(bg)%ndf%phi, &
                        phys%phase(bg)%uvwp%u,phys%phase(bg)%uvwp%v,phys%phase(bg)%uvwp%w,phys%phase(bg)%prop%rho)

      call calc_mip(eqn,prop,geom,.false.,dt)

    endif

    end associate
  enddo
    ! calculate mixture momentum in preparation for pressure solution


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

    associate(eqn=>phys%uvwp,prop=>phys%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

      call calc_coef_p(eqn,ap,anb,b,geom,prop)

      call set_a_0(phic,geom%ne+geom%nbf)
      call solve_gs('pc',phic,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

      pref=phic(1)! if closed system
      call adjust_pc(phic,pref,geom)
      call calc_grad(phic,eqn%gpc,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
      call update_uvwp(eqn,prop,geom,phic,geom%ne,geom%nf,geom%nbf)

    end associate

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
    associate(eqn=>phys%phase(n)%energy,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

    call calc_grad(eqn%t,eqn%gt,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)
    call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_energy(eqn,ap,anb,b,geom,prop,dt)

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

    n=phys%iphase
    ! solve phase continuity equations
    do n=1,NPHASES
      phys%iphase=n
      if(n==GAS) then
        associate(eqn=>phys%phase(n)%vfr,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

        call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

        call calc_coef_vfr(eqn,ap,anb,b,geom,prop,dt,phys%phase(n)%ndf%phi,phys%phase(n)%energy%t,phys%phase(n)%prop%rho)

        call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

        end associate
      elseif(n==LIQUID) then

        call calc_algebraic_vfr(phys%phase,n,geom%ne)

      endif

      call calc_phase_mfr(phys%phase,phys%phase(n)%mfr,geom)
    enddo

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

  subroutine calc_phase_mfr(phase,mfr,geom)
   type(phase_t) :: phase(:)
   real :: mfr(*)
   type(geometry_t) :: geom
   integer :: p,e
   real :: mfr0

   do e=1,geom%ne
     mfr0=0.
     do p=1,NPHASES
       mfr0=mfr(e)+phase(p)%prop%rho(e)*phase(p)%vfr%phi(e)
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

  subroutine solve_mfr(phys,geom) !eqn,prop,geom,dt,nit,ap,anb,b,phic)
    use mod_properties
    use mod_solver
    implicit none

    type(phys_t) :: phys
    type(geometry_t) :: geom
    real :: pref
    integer :: bg,p,c

    p=phys%iphase
    do c=1,NCOMPONENTS
      if(.not. phys%phase(p)%COMPONENT_SET(c)) cycle
      associate(eqn=>phys%phase(p)%component(c)%mfr,prop=>phys%phase(p)%component(c)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)
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
    integer :: bg,n

    n=phys%iphase
    associate(eqn=>phys%scalar,prop=>phys%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

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

  do n=1,NPHASES
    phys%iphase=n
    associate(eqn=>phys%phase(n)%ndf,prop=>phys%phase(n)%prop,ap=>phys%ap,anb=>phys%anb,b=>phys%b,phic=>phys%phic,nit=>phys%nit,dt=>phys%dt)

    !call calc_grad(eqn%phi,eqn%grad,geom%xc,geom%yc,geom%zc,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf)

    call calc_coef_ndf(eqn,ap,anb,b,geom,prop,dt)

    call solve_gs(eqn%name,eqn%phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)

    end associate
  enddo

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

    ne=geom%ne
  do n=1,NPHASES
    phys%iphase=n
    associate(prop=>phys%phase(n)%prop)

    ! calculate density
    if(iphase==GAS) then
      call calc_pvap(prop%pvap,prop%pvap_coef,phys%phase(n)%energy%t,ne)
      call calc_rho_idea(prop%rho,prop%rgas,phys%uvwp%p,phys%phase(n)%energy%t,ne)
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
