module mod_solver
  use mod_util
  use mod_subdomains
  implicit none

    character(len=*), parameter :: oformat = "(5x,A16,x,i5,x,15x,es9.3e2,3x,es9.3e2,3x,es9.3e2)"
 contains
  pure function matinv3(A) result(B)
  implicit none
  real, intent(in) :: A(3,3)   !! Matrix
  real             :: B(3,3)   !! Inverse matrix
  real             :: det,detinv

  det = (  A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  if(abs(det)>tiny(det)) then
    detinv=1._8/det

    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  else
    B=0._8
    det=A(1,1)+A(2,2)+A(3,3)
    detinv=1._8/det
    B(1,1)=detinv
    B(2,2)=detinv
    B(3,3)=detinv
  endif
  end function

  subroutine calc_grad(phi,grad,xc,yc,zc,ef2nb_idx,ef2nb,ne,nf,nbf)
    real, dimension(3,ne+nbf) :: grad
    real, dimension(ne+nbf) :: phi
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf
    integer :: e,idx,enb,lfnb
    real :: ga11,ga12,ga13,ga22,ga23,ga33,A(3,3),B(3,1),C(3,1),A_inv(3,3),dphi,dr(3),wt,xc(ne+nbf),yc(ne+nbf),zc(ne+nbf)

    do e=1,ne
      A=0.
      grad(:,e)=0.
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        dr = [xc(enb)-xc(e),yc(enb)-yc(e),zc(enb)-zc(e)]

        wt = 1./sum(dr**2)

        dphi=phi(enb)-phi(e)
        grad(:,e)=grad(:,e)+wt*dphi*dr

        A(1,1)=A(1,1)+wt*dr(1)*dr(1)
        A(1,2)=A(1,2)+wt*dr(1)*dr(2)
        A(1,3)=A(1,3)+wt*dr(1)*dr(3)
        A(2,2)=A(2,2)+wt*dr(2)*dr(2)
        A(2,3)=A(2,3)+wt*dr(2)*dr(3)
        A(3,3)=A(3,3)+wt*dr(3)*dr(3)

      end do

      A(2,1)=A(1,2)
      A(3,1)=A(1,3)
      A(3,2)=A(2,3)

      A_inv = matinv3(A)

      B(:,1)=grad(:,e)
      C=matmul(A_inv,B)
      grad(:,e)=C(:,1)

    enddo

  end subroutine

  subroutine internal_extrapolation_grad(phi,e,grad,xc,yc,zc,ef2nb_idx,ef2nb,ne,nf,nbf)
    real, dimension(3) :: grad
    real, dimension(ne+nbf) :: phi
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf
    integer :: e,idx,enb,lfnb
    real :: ga11,ga12,ga13,ga22,ga23,ga33,A(3,3),B(3,1),C(3,1),A_inv(3,3),dphi,dr(3),wt,xc(ne+nbf),yc(ne+nbf),zc(ne+nbf)

      grad=0.
      A=0.
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        if(lfnb==0) cycle
        dr = [xc(enb)-xc(e),yc(enb)-yc(e),zc(enb)-zc(e)]

        wt = 1./sum(dr**2)

        dphi=phi(enb)-phi(e)
        grad=grad+wt*dphi*dr

        A(1,1)=A(1,1)+wt*dr(1)*dr(1)
        A(1,2)=A(1,2)+wt*dr(1)*dr(2)
        A(1,3)=A(1,3)+wt*dr(1)*dr(3)
        A(2,2)=A(2,2)+wt*dr(2)*dr(2)
        A(2,3)=A(2,3)+wt*dr(2)*dr(3)
        A(3,3)=A(3,3)+wt*dr(3)*dr(3)

      end do

      A(2,1)=A(1,2)
      A(3,1)=A(1,3)
      A(3,2)=A(2,3)

      A_inv = matinv3(A)

      B(:,1)=grad
      C=matmul(A_inv,B)
      grad=C(:,1)

  end subroutine

    subroutine multi_subdomain_solver(cname,subdomain,intf,geom,ap,anb,b,phi,nsubd,nit,res_i_tot,res_f_tot,res_max_tot)
    use mod_solver_cuda
    use mod_util_cuda
    use omp_lib
    use cudafor
    implicit none
    integer :: nit,nsubd
    character(len=*) :: cname
    real , dimension(*):: ap,anb,b,phi
    type(geometry_t) :: geom
    type(subdomain_t) :: subdomain(nsubd)
    type(intf_t) :: intf(nsubd,nsubd)
    integer :: c ,i,j,it,idxout,devid
    character(len=32) :: name
    real :: res,res_max,res_i_tot,res_max_tot,res_f_tot,res_target,istat,starttime,endtime
    !real :: start,finish
    c = omp_get_thread_num()+1

     !assemble subdomain coef
   ! if (subdomain(c)%prefetch==1) then
  !  istat = cudaMemPrefetchAsync(subdomain(c)%phic, subdomain(c)%ne+subdomain(c)%nbf, cudaCpuDeviceId, subdomain(c)%stream)
  !  end if

!$omp barrier
!$omp single
    call assemble_coef(subdomain,geom,ap,anb,b,phi,nsubd)
!$omp endsingle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! \to be added. call cpu_time()
!$omp master
call cpu_time(starttime)
!$omp end master
!!$omp do
  !  do c=1,nsubd
  call initiate_device_data(subdomain(c))
 !     if(subdomain(c)%ARCH==GPU)
  !  enddo
!!$omp enddo
  !   solve
    res_i_tot=0.
    res_f_tot=0.
    res_max_tot=0.
    ! calculate initial residual

!!$omp do
!    do c=1,nsubd
      if(subdomain(c)%ARCH==CPU) then
        call calc_residual(subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
                         subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,res,res_max)

      elseif(subdomain(c)%ARCH==GPU) then
        call calc_residual_gpu(subdomain(c),res,res_max)
        !call calc_residual(subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
         !               subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,res,res_max)
      end if

!$omp critical
      res_i_tot=res_i_tot+res**2
      res_max_tot=max(res_max_tot,res_max)
!$omp end critical
 !   end do
!!$omp enddo
!$omp barrier
!$omp single
    !print *, subdomain(1)%b(5)
    res_i_tot=sqrt(res_i_tot/nsubd)
    ! loop until convergence criteria is fulfilled
    res_f_tot=res_i_tot
!     print * , res_i_tot
!      print *,  res_max_tot
!$omp end single
!!$omp master
!            call cpu_time(start)
!!$omp end master
    it = 0
    res_target=res_i_tot/10.! reduce residual 1 order of magnitude
    if(cname(1:2)=='pc') res_target=res_i_tot/10.
    do while(it<nit)
     !   print *, it
  !    do while(it<nit)! .and. res_f_tot>res_target)

      ! smooth each subdomain solution 2 iterations

!    !$omp single
!    print * , 'Phi(5) host',subdomain(1)%phic(5)
!    !$omp end single

!!$omp do
 !     do c=1,nsubd
        if(subdomain(c)%ARCH==CPU) then
          call smoother_gs(cname,subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
                         subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,2)
        elseif(subdomain(c)%ARCH==GPU) then

            call smoother_jacobi_gpu(subdomain(c))

            !if (subdomain(c)%prefetch==1) then
              !  call copy_d2h(subdomain(c))
           ! end if
        end if
 !     end do
!!$omp enddo

 !$omp barrier

      ! update interface halos
!!$omp do
   ! do i=1,nsubd
!
       i =c
      ! j=3
      ! print *,i
        do j=i+1,nsubd
        if (intf(i,j)%ncs==0) cycle
       ! print *, 'updated'
      ! call update_halos(intf(i,j),subdomain(i),subdomain(j),geom)
       call update_halos_gpu(intf(i,j),subdomain(i),subdomain(j),geom) ! c is thread Id


       ! end do
            !subdomain(c)%phic0 = subdomain(c)%phic
      end do
!!$omp enddo

 !$omp barrier
 ! if (subdomain(c)%prefetch==1) then
                !call copy_h2d(subdomain(c))
    !        endif
!$omp barrier


      ! calculate residuals every 10 iteration
      if(mod(it,10)==0) then
!!$omp do
!       do c=1,nsubd
          if(subdomain(c)%ARCH==CPU) then
            call calc_residual(subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
                         subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,res,res_max)

          elseif(subdomain(c)%ARCH==GPU) then
            call calc_residual_gpu(subdomain(c),res,res_max)
            !call calc_residual(subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
             !            subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,res,res_max)

          end if


!$omp critical
          res_f_tot=res_f_tot+res**2
          res_max_tot=max(res_max_tot,res_max)
!$omp end critical
!        end do
!!$omp enddo
!$omp barrier
!$omp single
        res_f_tot=sqrt(res_f_tot/nsubd)
!$omp end single
      endif

      it=it+1


    end do
    !$omp barrier
           !   if (subdomain(c)%prefetch==1) then
           !     call copy_d2h(subdomain(c))
           ! endif
!$omp barrier
!$omp master
call cpu_time(endtime)
geom%telap = geom%telap+endtime-starttime
!$omp end master

!$omp single

    name=trim(cname)
    write(*,oformat) name,it,res_i_tot,res_f_tot,res_max_tot
    ! update main domain phi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_phi(subdomain,phi,geom,nsubd)
!$omp endsingle
  end subroutine


!  subroutine multi_subdomain_solver(cname,subdomain,intf,geom,ap,anb,b,phi,nsubd,nit,res_i_tot,res_f_tot,res_max_tot)
!    use mod_solver_cuda
!    use mod_util_cuda
!    use omp_lib
!    implicit none
!    integer :: nit,nsubd
!    character(len=*) :: cname
!    real , dimension(*):: ap,anb,b,phi
!    type(geometry_t) :: geom
!    type(subdomain_t) :: subdomain(nsubd)
!    type(intf_t) :: intf(nsubd,nsubd)
!    integer :: c ,i,j,it,idxout
!    character(len=32) :: name
!    real :: res,res_max,res_i_tot,res_max_tot,res_f_tot,res_target,istat
!    real :: start,finish
!    integer :: tx,gx,sz,length
!    type(dim3) :: grid,block
!
!
!
!    c = omp_get_thread_num()+1
!
!
!
!!$omp single
!    call assemble_coef(subdomain,geom,ap,anb,b,phi,nsubd)
!!$omp endsingle
!
!!call allocate_device_data(subomain(c))
!
!!print *, subdomain(c)%stream
!!!$omp barrier
!!    istat = cudaMemPrefetchAsync(subdomain(c)%ap, subdomain(c)%ne, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%anb, 2*subdomain(c)%nf-subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%b, subdomain(c)%ne, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%phic, subdomain(c)%ne+subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%phic0, subdomain(c)%ne+subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%ef2nb, 2*subdomain(c)%nf-subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!    istat = cudaMemPrefetchAsync(subdomain(c)%ef2nb_idx, subdomain(c)%ne+1, subdomain(c)%device_num, subdomain(c)%stream)
!
!    ! solve
!    res_i_tot=0.
!    res_f_tot=0.
!    res_max_tot=0.
!    ! calculate initial residual
!
!    !call calc_residual_gpu(subdomain(c),res,res_max)
!
!    tx = 512
!    gx = subdomain(c)%ne/tx+1
!    grid = dim3(gx,1,1)
!    block = dim3(tx,1,1)
!
!
!    call calc_residual_kernel<<<grid,block,0,subdomain(c)%stream>>>(subdomain(c)%phic,&
!                        subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,&
!                        subdomain(c)%ef2nb,subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,subdomain(c)%r)!,subdomain%res_max)
!
!    !print *, subdomain%stream
!    !istat = cudaStreamSynchronize(subdomain%stream)
!    tx=512*2
!    grid = dim3(1,1,1)
!    block = dim3(tx,1,1)
!    ! size may potentially cause trouble as it is dependednt on thread num because of limit on shared memory
!    sz=tx*sizeof(res)*2! number of bytes for shared memory
!
!     ! call residual_reduction_unified<<<grid,block,sz,subdomain%stream>>>(subdomain%r,subdomain%ne,tx,subdomain%res)
!    call residual_reduction<<<grid,block,sz,subdomain(c)%stream>>>(subdomain(c)%r,subdomain(c)%ne,tx,subdomain(c)%res(1),subdomain(c)%res_max(1))
!
!   ! istat = cudaStreamSynchronize(subdomain%stream)
!    !istat= cudaDeviceSynchronize()
!   ! istat = cudaMemPrefetchAsync(subdomain(c)%res, 2, cudaCpuDeviceId, subdomain(c)%stream)
!
!    res = subdomain(c)%res(1)!(1)
!    res_max = subdomain(c)%res_max(1) !(2)
!
!
!!$omp critical
!      res_i_tot=res_i_tot+res**2
!      res_max_tot=max(res_max_tot,res_max)
!!$omp end critical
!!$omp barrier
!!$omp single
!    res_i_tot=sqrt(res_i_tot/nsubd)
!    res_f_tot=res_i_tot
!!$omp end single
!
!    it = 0
!    res_target=res_i_tot/10.! reduce residual 1 order of magnitude
!    if(cname(1:2)=='pc') res_target=res_i_tot/10.
!    do while(it<100)
!
!    !    call smoother_jacobi_gpu(subdomain(c))
!
!    tx = 256
!    gx = subdomain(c)%ne/tx+1
!    grid = dim3(gx,1,1)
!    block = dim3(tx,1,1)
!
!
!    call solve_jacobi_kernel<<<grid,block,0,subdomain(c)%stream>>>(subdomain(c)%phic,subdomain(c)%phic0,&
!                        subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,&
!                        subdomain(c)%ef2nb,subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf)
!
!  !  istat = cudaMemPrefetchAsync(subdomain(c)%phic, subdomain(c)%ne+subdomain(c)%nbf, cudaCpuDeviceId, subdomain(c)%stream)
!
!
! !$omp barrier
!! !$omp master
!!            call cpu_time(start)
!!!$omp end master
!        i=c
!        do j=i+1,nsubd
!          call update_halos(intf(i,j),subdomain(i),subdomain(j),geom)
!        end do
!
!!  istat = cudaMemPrefetchAsync(subdomain(c)%phic, subdomain(c)%ne+subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!  istat = cudaMemPrefetchAsync(subdomain(c)%phic0, subdomain(c)%ne+subdomain(c)%nbf, subdomain(c)%device_num, subdomain(c)%stream)
!!
!
!      !end do
!! ! calculate residuals every 10 iteration
!      if(mod(it,10)==0) then
!!!$omp do
!!       do c=1,nsubd
!
!            call calc_residual_gpu(subdomain(c),res,res_max)
!            !call calc_residual(subdomain(c)%phic,subdomain(c)%ap,subdomain(c)%anb,subdomain(c)%b,subdomain(c)%ef2nb_idx,subdomain(c)%ef2nb,&
!             !            subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf,res,res_max)
!!$omp critical
!          res_f_tot=res_f_tot+res**2
!          res_max_tot=max(res_max_tot,res_max)
!
!!$omp end critical
!!        end do
!!!$omp enddo
!!$omp barrier
!!$omp single
!        res_f_tot=sqrt(res_f_tot/nsubd)
!!$omp end single
!      endif
!
!
!
!it=it+1
!    end do
!
!
!!$omp single
!    name=trim(cname)
!    write(*,oformat) name,it,res_i_tot,res_f_tot,res_max_tot
!    ! update main domain phi
!
!    call update_phi(subdomain,phi,geom,nsubd)
!!$omp endsingle
!  end subroutine

subroutine smoother_jacobi_cpu_debug(cname,phi,phi0,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,nit)
    implicit none
    character(len=*) :: cname
    real, dimension(ne+nbf) :: phi,phi0
    real, dimension(ne) :: ap,b
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf,nit,it
    integer :: e,idx,enb,lfnb
    real :: sumnb,sor



    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi0(enb)
      end do
      phi(e)=(sumnb)/ap(e)
      !phi0(e)=phi(e)
    end do
end subroutine


  subroutine smoother_gs(cname,phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,nit)
    implicit none
    character(len=*) :: cname
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf,nit,it
    integer :: e,idx,enb,lfnb
    real :: sumnb,sor

    sor=1.
    if(cname(1:2)=='pc') sor=1.02

  do it=1,nit
    ! forward sweep
    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      phi(e)=(sumnb+(sor-1.)*ap(e)*phi(e))/ap(e)/sor
    end do

    ! backward sweep
    do e=ne,1,-1
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      phi(e)=(sumnb+(sor-1.)*ap(e)*phi(e))/ap(e)/sor
    end do

  enddo

  end subroutine

  subroutine calc_residual(phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,res,res_max)
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf,nit,it
    integer :: e,idx,enb,lfnb
    real :: sumnb,res,res_max,r

    res=0.
    res_max=0.
    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      r=sumnb-ap(e)*phi(e)
      res_max=max(r,res_max)
      res=res+r**2
    end do

    res=sqrt(res/ne)

  end subroutine

  subroutine solve_gs(cname,phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,nit)
    implicit none
    character(len=*) :: cname
    character(len=16) :: name
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: ne,nf,nbf,nit,it
    integer :: e,idx,enb,lfnb
    real :: sumnb,res_i,res_f,res_target,res_max,r
    real :: sor

    sor=1.
    if(cname(1:2)=='pc') sor=1.02

    res_i=0.! initial residual
    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      r=sumnb-ap(e)*phi(e)
      res_i=res_i+r**2
    end do
    res_i=sqrt(res_i/ne)

  it=0
  res_f=res_i
  res_target=res_i/10.! reduce residual 1 order of magnitude
  if(cname(1:2)=='pc') res_target=res_i/10.
  do while(it<nit .and. res_f>res_target)
    it=it+1
    ! forward sweep
    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      phi(e)=(sumnb+(sor-1.)*ap(e)*phi(e))/ap(e)/sor
    end do

    ! backward sweep
    do e=ne,1,-1
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      phi(e)=(sumnb+(sor-1.)*ap(e)*phi(e))/ap(e)/sor
    end do

    res_max=0.
    res_f=0.! final residual
    do e=1,ne
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      r=abs(sumnb-ap(e)*phi(e))
      res_max=max(r,res_max)
      res_f=res_f+r**2
    end do
    res_f=sqrt(res_f/ne)
  enddo

    name=trim(cname)
    write(*,oformat) name,it,res_i,res_f,res_max

  end subroutine

  subroutine solve(cname,subdomain,intf,geom,ap,anb,b,phi,nsubd,nit)
    use omp_lib
    use mod_solver_cuda
    implicit none
    integer :: nit,nsubd,istat,tid
    integer(kind=cuda_stream_kind) :: mystream
    real :: startt,endt

    character(len=*) :: cname
    real , dimension(*):: ap,anb,b,phi
    type(geometry_t) :: geom
    type(subdomain_t) :: subdomain(nsubd)
    type(intf_t) :: intf(nsubd,nsubd)

    real :: res_i,res_f,res_max

!    if(nsubd==1) then
!      call solve_gs(cname,phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
!    else
    if(nsubd==1) then
    call cpu_time(startt)

    !  call solve_gs(cname,phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)



!   call solve_jacobi_db_gpu(cname,geom,phi,ap,anb,b,geom%ef2nb_idx,geom%ef2nb,geom%ne,geom%nf,geom%nbf,nit)
    !!!!!!!!! SINGLE GPU THING!!!


  call cpu_time(endt)
    geom%telap = geom%telap + (endt-startt)
    else

!$omp parallel num_threads(nsubd) shared(cname,subdomain, intf,geom,ap,anb,b,phi,nsubd,nit,res_i,res_f,res_max) default(private)
        tid = omp_get_thread_num()
        !sub_for_tid


        !istat = cudaStreamCreate(subdomain(tid+1)%stream)
        istat = cudaforSetDefaultstream(subdomain(tid+1)%stream)
       ! print *, 'first', tid,subdomain(tid+1)%stream
        call multi_subdomain_solver(cname,subdomain,intf,geom,ap,anb,b,phi,nsubd,nit,res_i,res_f,res_max)
!$omp endparallel
    end if
    !call cpu_time(endt)
   ! geom%telap = geom%telap+(endt-startt)
  end subroutine
end module
