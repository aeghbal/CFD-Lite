module mod_solver_cuda
  implicit none
  contains

  subroutine allocate_device_data(subdomain)
    use mod_util
    use mod_subdomains
    implicit none


    type(subdomain_t) :: subdomain
    integer :: istat,length,devid

  !  do devid=1,4
    istat = cudaSetDevice(subdomain%device_num)
!    length=size(subdomain%ef2nb)
!    allocate(subdomain%ef2nb_d(length))
!    istat = cudaMemcpyAsync(subdomain%ef2nb_d,subdomain%ef2nb,length,subdomain%stream)
!    length=size(subdomain%ef2nb_idx)
!    allocate(subdomain%ef2nb_idx_d(length)) !Updated this line. Previously - allocate(subdomain%ef2nb_d(length))
!    istat = cudaMemcpyAsync(subdomain%ef2nb_idx_d,subdomain%ef2nb_idx,length,subdomain%stream)
!  !  end do
!    subdomain%ef2nb_d = subdomain%ef2nb
!    subdomain%ef2nb_idx = subdomain%ef2nb_idx


    length=size(subdomain%ap)
    allocate(subdomain%ap_d(length))
    allocate(subdomain%b_d(length))
 !   allocate(subdomain%r_d(length))
!    allocate(subdomain%res_d)
!    allocate(subdomain%res_max_d)
    length=size(subdomain%anb)
    allocate(subdomain%anb_d(length))
!    length=size(subdomain%phic)
!    allocate(subdomain%phic_d(length))
!    allocate(subdomain%phic0_d(length))
    istat = cudaDeviceSynchronize()

   end subroutine

  subroutine initiate_device_data(subdomain)
    use mod_util
    use mod_subdomains
    implicit none

    type(subdomain_t) :: subdomain
    integer :: istat,devid

   ! print *, 'preintitated'
    istat = cudaSetDevice(subdomain%device_num)
    istat = cudaMemcpyAsync(subdomain%ap_d,subdomain%ap,size(subdomain%ap),subdomain%stream)
    istat = cudaMemcpyAsync(subdomain%anb_d,subdomain%anb,size(subdomain%anb),subdomain%stream)
    istat = cudaMemcpyAsync(subdomain%b_d,subdomain%b,size(subdomain%b),subdomain%stream)

    istat = cudaDeviceSynchronize()
   ! print *, 'intitated'

  end subroutine



  subroutine copy_d2h(subdomain)
    use mod_util
    use mod_subdomains
    use cudafor
    implicit none


    type(subdomain_t) :: subdomain
    integer :: istat,length

    istat = cudaSetDevice(subdomain%device_num)
    istat = cudaMemPrefetchAsync(subdomain%phic, subdomain%ne+subdomain%nbf, cudaCpuDeviceId, subdomain%stream)
 !   istat = cudaMemPrefetchAsync(subdomain%phic0, subdomain%ne+subdomain%nbf, cudaCpuDeviceId, subdomain%stream)

    istat = cudaDeviceSynchronize()

  end subroutine

  subroutine copy_h2d(subdomain)
    use mod_util
    use mod_subdomains
    implicit none


    type(subdomain_t) :: subdomain
    integer :: istat,length
!UM = 35
!UM +Prefetch = 42
!Um + Advice =25

    istat = cudaSetDevice(subdomain%device_num)
    istat = cudaMemPrefetchAsync(subdomain%phic, subdomain%ne+subdomain%nbf, subdomain%device_num, subdomain%stream)
   ! istat = cudaMemPrefetchAsync(subdomain%phic0, subdomain%ne+subdomain%nbf, subdomain%device_num, subdomain%stream)
    istat = cudaDeviceSynchronize()

  end subroutine

  subroutine calc_residual_gpu(subdomain,res,res_max)
    use mod_util
    use mod_subdomains
    use mod_util_cuda
    use omp_lib
    implicit none
    type(subdomain_t) :: subdomain
    integer :: istat,tx,gx,sz
    type(dim3) :: grid,block
    real :: res,res_max

    istat = cudaSetDevice(subdomain%device_num)

    tx = 512
    gx = subdomain%ne/tx+1
    grid = dim3(gx,1,1)
    block = dim3(tx,1,1)

    if (subdomain%prefetch==1) then
    istat = cudaMemPrefetchAsync(subdomain%phic, subdomain%ne+subdomain%nbf, subdomain%device_num, subdomain%stream)
    istat = cudaMemPrefetchAsync(subdomain%res, 2, subdomain%device_num, subdomain%stream)
    istat = cudaMemPrefetchAsync(subdomain%r, subdomain%ne, subdomain%device_num, subdomain%stream)
    end if
   ! istat = cudaDeviceSynchronize()

!    call calc_residual_kernel<<<grid,block,0,subdomain%stream>>>(subdomain%phic,&
!                        subdomain%ap,subdomain%anb,subdomain%b,subdomain%ef2nb_idx,&
!                        subdomain%ef2nb,subdomain%ne,subdomain%nf,subdomain%nbf,subdomain%r)!,subdomain%res_max)

    call calc_residual_kernel<<<grid,block,0,subdomain%stream>>>(subdomain%phic,&
                        subdomain%ap_d,subdomain%anb_d,subdomain%b_d,subdomain%ef2nb_idx,&
                        subdomain%ef2nb,subdomain%ne,subdomain%nf,subdomain%nbf,subdomain%r)


    !print *, subdomain%stream
    !istat = cudaStreamSynchronize(subdomain%stream)
    tx=512*2
    grid = dim3(1,1,1)
    block = dim3(tx,1,1)
    ! size may potentially cause trouble as it is dependednt on thread num because of limit on shared memory
    sz=tx*sizeof(res)*2! number of bytes for shared memory

     ! call residual_reduction_unified<<<grid,block,sz,subdomain%stream>>>(subdomain%r,subdomain%ne,tx,subdomain%res)
    call residual_reduction<<<grid,block,sz,subdomain%stream>>>(subdomain%r,subdomain%ne,tx,subdomain%res(1),subdomain%res_max(1))

   ! istat = cudaStreamSynchronize(subdomain%stream)
    istat= cudaDeviceSynchronize()
    !istat = cudaSetDevice(subdomain%device_num)
    istat = cudaMemPrefetchAsync(subdomain%res, 1, cudaCpuDeviceId, subdomain%stream)
    istat = cudaMemPrefetchAsync(subdomain%res_max, 1, cudaCpuDeviceId, subdomain%stream)


    res = subdomain%res(1)!(1)
    res_max = subdomain%res_max(1) !(2)


  end subroutine

  subroutine smoother_jacobi_gpu(subdomain)
    use mod_util
    use mod_subdomains
    use mod_util_cuda
    implicit none
    type(subdomain_t) :: subdomain
    integer :: istat,tx,gx,sz,length
    type(dim3) :: grid,block
    real :: res,res_max

    istat = cudaSetDevice(subdomain%device_num)

    if (subdomain%prefetch==1) then
    istat = cudaMemPrefetchAsync(subdomain%phic, subdomain%ne+subdomain%nbf, subdomain%device_num, subdomain%stream)
    istat = cudaMemPrefetchAsync(subdomain%phic0, subdomain%ne+subdomain%nbf, subdomain%device_num, subdomain%stream)
    end if


    tx = 256
    gx = subdomain%ne/tx+1
    grid = dim3(gx,1,1)
    block = dim3(tx,1,1)

    call solve_jacobi_kernel<<<grid,block,0,subdomain%stream>>>(subdomain%phic,subdomain%phic0,&
                        subdomain%ap_d,subdomain%anb_d,subdomain%b_d,subdomain%ef2nb_idx,&
                        subdomain%ef2nb,subdomain%ne,subdomain%nf,subdomain%nbf)



    istat = cudaDeviceSynchronize()

  end subroutine

  subroutine update_halos_gpu(intf,subdomain1,subdomain2,geom)
    use mod_subdomains
    use cudafor
    use mod_util_cuda
    use omp_lib

    implicit none

    type(subdomain_t) :: subdomain1,subdomain2
    type(intf_t) :: intf
    type(geometry_t) :: geom
    type(dim3) :: grid,block
    integer :: istat,tx,gx,sz_id1,sz_id2,ddd,dddd!,g1(intf%ncs),g2(intf%ncs),gnb1(intf%ncs),gnb2(intf%ncs)


    integer :: nsubd,th_id!cs,g1,g2,lf1,lf2,gnb1,gnb2,tmp,idx,

    if(intf%c1==intf%c2) return
    tx = 256
    gx = intf%ncs/tx+1
    grid = dim3(gx,1,1)
    block = dim3(tx,1,1)
     !   print *, th_id,nsubd

    sz_id1=size(intf%index1)
    sz_id2=size(intf%index2)
    !print *, sz_id1,sz_id2,intf%ncs

     istat = cudaSetDevice(subdomain1%device_num)


!    call update_halos_kernel<<<grid,block,0,subdomain1%stream>>>(subdomain1%ef2nb_idx_d,subdomain1%ef2nb_d, intf%g1_d,intf%gnb1_d,&
!                    intf%index1_d,intf%ncs,subdomain1%ne,subdomain1%nf,subdomain1%nbf,sz_id1)

  !  istat = cudaSetDevice(subdomain2%device_num)

!    call update_halos_kernel<<<grid,block,0,subdomain2%stream>>>(subdomain2%ef2nb_idx_d,subdomain2%ef2nb_d, intf%g2_d,intf%gnb2_d,&
!                    intf%index2_d,intf%ncs,subdomain2%ne,subdomain2%nf,subdomain2%nbf,sz_id2)

    call update_halos_kernel<<<grid,block,0,subdomain1%stream>>>(subdomain1%ef2nb_idx,subdomain1%ef2nb, subdomain1%phic,&
                        subdomain2%ef2nb_idx,subdomain2%ef2nb, subdomain2%phic,intf%index1_d,intf%index2_d,intf%ncs,&
                        subdomain1%ne,subdomain1%nf,subdomain1%nbf,subdomain2%ne,subdomain2%nf,subdomain2%nbf,sz_id1,sz_id2)
!
    istat = cudaDeviceSynchronize()


  end subroutine

!     subroutine solve_jacobi_db_gpu(cname,geom,phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,nit)
!    use cudafor
!    use mod_util
!    use mod_util_cuda
!    use mod_mg_lvl
!    implicit none
!    character(len=*) :: cname
!    character(len=16) :: name
!    type(geometry_t) :: geom
!
!    real, dimension(ne+nbf) :: phi,phi0
!    real, dimension(ne) :: ap,b
!    real, dimension(2*nf-nbf) :: anb
!    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
!    integer :: ne,nf,nbf,nit,it
!    integer :: e,idx,enb,lfnb
!    real :: res_i,res_f,res_target
!!    real,device, allocatable, dimension(:) :: ap_d,anb_d,b_d,phic_d,phic0_d,r_d
!!    real,allocatable, device :: res_d, res_max_d
!!    integer, device, allocatable :: ef2nb_d(:),ef2nb_idx_d(:)
!    type(dim3) :: grid,block
!    integer :: istat,tx,gx,sz,devvv
!    real :: res,res_max
!    character(len=*), parameter :: oformat = "(5x,A16,x,i5,x,15x,es9.3e2,3x,es9.3e2,3x,es9.3e2)"
!
!
!
!    phi0=phi
!    res_i=0.! initial residual
!
!!    geom%ap=ap
!!    geom%b=b
!!    geom%anb=anb
!!    geom%phic=phic
!
!    tx = 512
!    gx = ne/tx+1
!    grid = dim3(gx,1,1)
!    block = dim3(tx,1,1)
!
!
!     !  print *, 'com['
!
!   ! print *,'in kernel', omp_get_thread_num() , subdomain%stream
!
!    call calc_residual_kernel<<<grid,block,0,0>>>(geom%phic, geom%ap,geom%anb,geom%b,geom%ef2nb_idx_d,geom%ef2nb_d,ne,nf,nbf,geom%r_d)!,subdomain%res_max)
!    istat = cudaDeviceSynchronize()
!
!  !  print *,'inres', subdomain%stream
!    tx=512*2
!    grid = dim3(1,1,1)
!    block = dim3(tx,1,1)
!    ! size may potentially cause trouble as it is dependednt on thread num because of limit on shared memory
!    sz=tx*sizeof(res)*2! number of bytes for shared memory
!
!    call residual_reduction<<<grid,block,sz,0>>>(geom%r_d,ne,tx,geom%res_d(1),geom%res_max_d(1))
!    istat = cudaDeviceSynchronize()
!
!    istat = cudaMemPrefetchAsync(geom%res_d, 1, cudaCpuDeviceId, 0)
!    istat = cudaMemPrefetchAsync(geom%res_max_d, 1, cudaCpuDeviceId, 0)
!    res = geom%res_d(1)
!    res_max = geom%res_max_d(1)
!
!    !!!!!!!!!!!!!!!!!
!    res_i=sqrt(res/ne)
!    res_f = res_i
!
!  it=0
!  res_f=res_i
!  res_target=res_i/10.! reduce residual 1 order of magnitude
!  if(cname(1:2)=='pc') res_target=res_i/10.
!  do while(it<nit)! .and. res_f>res_target)
!    it=it+1
!
!     tx = 256
!    gx = ne/tx+1
!    grid = dim3(gx,1,1)
!    block = dim3(tx,1,1)
!
!   ! print *,'in jacobi', subdomain%stream
!    call solve_jacobi_kernel<<<grid,block,0,0>>>(geom%phic,geom%phic0_d,geom%ap,geom%anb,geom%b,geom%ef2nb_idx_d,geom%ef2nb_d,ne,nf,nbf)
!
!
!    istat = cudaDeviceSynchronize()
!
!    res_max=0.
!    res_f=0.! final residual
!
!    tx = 512
!    gx = ne/tx+1
!    grid = dim3(gx,1,1)
!    block = dim3(tx,1,1)
!
!
!    !   print *, 'com[',devvv, subdomain%device_num
!
!   ! print *,'in kernel', omp_get_thread_num() , subdomain%stream
!    if(mod(it,10)==0) then
!
!    call calc_residual_kernel<<<grid,block,0,0>>>(geom%phic, geom%ap,geom%anb,geom%b,geom%ef2nb_idx_d,geom%ef2nb_d,ne,nf,nbf,geom%r_d)!,subdomain%res_max)
!    istat = cudaDeviceSynchronize()
!
!  !  print *,'inres', subdomain%stream
!    tx=512*2
!    grid = dim3(1,1,1)
!    block = dim3(tx,1,1)
!    ! size may potentially cause trouble as it is dependednt on thread num because of limit on shared memory
!    sz=tx*sizeof(res)*2! number of bytes for shared memory
!
!    call residual_reduction<<<grid,block,sz,0>>>(geom%r_d,ne,tx,geom%res_d(1),geom%res_max_d(1))
!    istat = cudaDeviceSynchronize()
!
!    istat = cudaMemPrefetchAsync(geom%res_d, 1, cudaCpuDeviceId, 0)
!    istat = cudaMemPrefetchAsync(geom%res_max_d, 1, cudaCpuDeviceId, 0)
!    res = geom%res_d(1)
!    res_max = geom%res_max_d(1)
!    !!!!!!!!!!!!!!!!!
!    res_f=sqrt(res/ne)
!    endif
!
!  enddo
!
!    name=trim(cname)
!    write(*,oformat) name,it,res_i,res_f,res_max
!
!  end subroutine

end module
