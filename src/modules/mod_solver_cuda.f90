module mod_solver_cuda
  implicit none

  contains

  subroutine allocate_device_data(subdomain)
    use mod_util
    use mod_subdomains
    implicit none

    type(subdomain_t) :: subdomain
    integer :: istat,length

    istat = cudaSetDevice(subdomain%device_num)

    length=size(subdomain%ef2nb)
    allocate(subdomain%ef2nb_d(length))
    istat = cudaMemcpyAsync(subdomain%ef2nb_d,subdomain%ef2nb,length,subdomain%stream)
    length=size(subdomain%ef2nb_idx)
    allocate(subdomain%ef2nb_d(length))
    istat = cudaMemcpyAsync(subdomain%ef2nb_idx_d,subdomain%ef2nb_idx,length,subdomain%stream)
    length=size(subdomain%ap)
    allocate(subdomain%ap_d(length))
    allocate(subdomain%b_d(length))
    allocate(subdomain%res_d(length))
    length=size(subdomain%anb)
    allocate(subdomain%anb_d(length))
    length=size(subdomain%phic)
    allocate(subdomain%phic_d(length))

  end subroutine

  subroutine initiate_device_data(subdomain)
    use mod_util
    use mod_subdomains
    implicit none

    type(subdomain_t) :: subdomain
    integer :: istat,length

    istat = cudaSetDevice(subdomain%device_num)
    istat = cudaMemcpyAsync(subdomain%ap_d,subdomain%ap,size(subdomain%ap),subdomain%stream)
    istat = cudaMemcpyAsync(subdomain%anb_d,subdomain%anb,size(subdomain%anb),subdomain%stream)
    istat = cudaMemcpyAsync(subdomain%b_d,subdomain%b,size(subdomain%b),subdomain%stream)
    istat = cudaMemcpyAsync(subdomain%phic_d,subdomain%phic,size(subdomain%phic),subdomain%stream)
  end subroutine

  subroutine calc_residual_gpu(subdomain,res,res_max)
    use mod_util
    use mod_subdomains
    use mod_util_cuda
    implicit none
    type(subdomain_t) :: subdomain
    integer :: istat,tx,gx,sz
    type(dim3) :: grid,block
    real :: res,res_max

    istat = cudaSetDevice(subdomain%device_num)

    tx = 256
    gx = cl%ne/tx+1
    grid = dim3(gx,1,1)
    block = dim3(tx,1,1)

    call calc_residual_kernel<<<grid,block,0,subdomain%stream>>>(subdomain(c)%phic_d,&
                        subdomain(c)%ap_d,subdomain(c)%anb_d,subdomain(c)%b_d,subdomain(c)%ef2nb_idx_d,&
                        subdomain(c)%ef2nb_d,subdomain(c)%ne,subdomain(c)%nf,subdomain(c)%nbf)

    tx=512
    grid = dim3(1,1,1)
    block = dim3(tx,1,1)
    sz=tx*sizeof(res)! number of bytes for shared memory

    call residual_reduction<<<grid,block,sz,subdomain%stream>>>(subdomain%res_d,subdomain%ne,tx,res,res_max)

  end subroutine


end module
