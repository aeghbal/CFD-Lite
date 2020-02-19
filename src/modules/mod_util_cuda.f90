module mod_util_cuda
  use cudafor
  implicit none

  contains

  subroutine calc_residual_kernel(phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,res)
    integer, value :: ne,nf,nbf
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b,res
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: e,idx,enb,lfnb
    real :: sumnb


    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
    do while(e<=ne)
      res(e)=0.
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      res(e)=sumnb-ap(e)*phi(e)
    end do

  end subroutine

  subroutine residual_reduction(r,length,nthrd,res,res_max)
    real, value :: res
    integer, value :: length,nthrd
    real, dimension(length) :: r
    real, shared, dimension(nthrd,2) :: sh! length equal to num threads X 8 bytes
    integer :: tid,bx
    ! residual reduction
    tid = threadIdx%x
    bx = blockDim%x

    e = tid
    sh(tid,:) = constRealType(0.0)

    do while(e <= length)
      sh(tid,1) = sh(tid,1) + abs(r(e))
      sh(tid,2) = max(sh(tid,2),abs(r(e)))
      e = e + bx
    enddo

    call syncthreads()

    bx = ishft(bx,-1)
    do while(bx > 0)
      if(tid <= bx) then
        sh(tid,1) = sh(tid,1) + sh(tid+bx,1)
        sh(tid,2) = max(sh(tid,2),sh(tid+bx,2))
      endif
      bx = ishft(bx,-1)
      call syncthreads()
    enddo

    if(tid == 1) red = sh(1)

  end subroutine
