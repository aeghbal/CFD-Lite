module mod_util_cuda
  use cudafor
  implicit none



  integer, parameter :: sec_lvl=3,num_face_bits=5
  contains
  attributes(device) subroutine get_idx_gpu(hs,lvl,g,s)
    implicit none
    integer :: hs,g,s,lvl

    s=iand(hs,2**(min(2*lvl+num_face_bits,16))-1)
    g=ishft(hs,-(min(2*lvl+num_face_bits,16)))
    !print *, s,g

  end subroutine


  attributes(global) subroutine calc_residual_kernel(phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,r)
    integer, value :: ne,nf,nbf
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b,r
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: e,idx,enb,lfnb,idxout,stride

    real :: sumnb

    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
    stride = blockDim%x * gridDim%x

 !   do while(e<=ne)
    Do idxout = e,ne,stride

      r(e)=0.
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx_gpu(ef2nb(idx,1),0,enb,lfnb)
        sumnb=sumnb+anb(idx)*phi(enb)
      end do
      r(e)=sumnb-ap(e)*phi(e)
      end do
     ! if (e==1)   print *, 'In residual kernel'
  end subroutine

  attributes(global) subroutine calc_residual_kernel_unified(phi,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf,r)
    integer, value :: ne,nf,nbf
    real, dimension(ne+nbf) :: phi
    real, dimension(ne) :: ap,b,r
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
    integer :: e,idx,enb,lfnb,idxout,stride

    real :: sumnb!,res(2)!,res_max(*)

    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
    stride = blockDim%x * gridDim%x

    Do idxout = e,ne,stride
      r(e)=0.
      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx_gpu(ef2nb(idx,1),0,enb,lfnb)
              sumnb=sumnb+anb(idx)*phi(enb)
      end do
      r(e)=sumnb-ap(e)*phi(e)
      end do
!    if (e==1) then
!    res(1)=0.
!    res(2)=0.!res_max=0.
!    do idxout=1,ne
!      res(1) = res(1)+r(idxout)**2
!      res(2)=max(r(idxout),res(2))
!      !res_max=max(r(idxout),res_max)
!    end do
!        res(1)=sqrt(res(1)/ne)
!    end if

  end subroutine

  attributes(global) subroutine residual_reduction(r,length,nthrd,res,res_max)
    real :: res,res_max !Here
    integer, value :: length,nthrd
    real, dimension(length) :: r
    real, shared, dimension(nthrd,2) :: sh! length equal to num threads X 8 bytes
    integer :: tid,bx,e
    ! residual reduction
    tid = threadIdx%x+(blockIdx%x-1)*blockDim%x
    bx = blockDim%x * gridDim%x
    e = tid
    sh(tid,:) = 0.0

    do while(e <= length)
      sh(tid,1) = sh(tid,1) + (r(e))**2
      sh(tid,2) = max(sh(tid,2),(r(e)))
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

    if(tid == 1) then
        res = sqrt(sh(1,1)/length)
        res_max = sh(1,2)
       ! print *, 'In residual reduction ernel'
    end if

  end subroutine

  attributes(global) subroutine residual_reduction_unified(r,length,nthrd,res)
    real :: res(2) !Here
    integer, value :: length,nthrd
    real, dimension(length) :: r
    real, shared, dimension(nthrd,2) :: sh! length equal to num threads X 8 bytes
    integer :: tid,bx,e
    ! residual reduction
    tid = threadIdx%x+(blockIdx%x-1)*blockDim%x
    bx = blockDim%x * gridDim%x
    e = tid
    sh(tid,:) = 0.0

    do while(e <= length)
      sh(tid,1) = sh(tid,1) + (r(e))**2
      sh(tid,2) = max(sh(tid,2),(r(e)))
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

    if(tid == 1) then
        res(1) = sqrt(sh(1,1)/length)
        res(2) = sh(1,2)
    end if

  end subroutine



  attributes(global) subroutine solve_jacobi_kernel(phi,phi0,ap,anb,b,ef2nb_idx,ef2nb,ne,nf,nbf)
    integer, value :: ne,nf,nbf
    real, dimension(ne+nbf) :: phi,phi0
    real, dimension(ne) :: ap,b
    real, dimension(2*nf-nbf) :: anb
    integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf)!,2)!
    integer :: e,idx,enb,lfnb,bx,idxout
    real :: sumnb


    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
    bx = blockDim%x * gridDim%x

!    do while(e<=ne)
    Do idxout = e,ne,bx

      sumnb=b(e)
      do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
        call get_idx_gpu(ef2nb(idx),0,enb,lfnb)
!       call get_idx_gpu(ef2nb(idx,1),0,enb,lfnb)

        sumnb=sumnb+anb(idx)*phi0(enb)
      end do
      phi(e)=sumnb/ap(e)
   !   phi0(e)= phi(e) ! NOTE THIS Line didnt exist in earlier versions.


!if coefficient are positive sum of a and b is on right side for this and exn.
! multiple blocks test - by end of next week. multiple-gpu
! perf. unified. GPU code working in concert with OMP-CUDA - unified vs seperate
! BC
!sliding lid - unified vs
! mpi comm propagator mgr. robust comm. env.
    end do
  !  if (e==1)   print *, 'In Jacobi kernel'

  end subroutine

!  attributes(global) subroutine update_halos_kernel(ef2nb_idx1,ef2nb1, g1,gnb1,index1,ncs,ne1,nf1,nbf1,sz_id1)
!
!    integer :: idxout,lf1,tmp,idx!,lf2
!    integer :: e,bx
!    integer, value :: ne1,nf1,nbf1,ncs,sz_id1!,sz_id2,ne2,nf2,nbf2
!    integer :: index1(sz_id1),g1(ncs),gnb1(ncs)!,gnb2(ncs),g2(ncs),index2(sz_id2)
!    !real :: phic1(ne1+nbf1), phic2(ne2+nbf2)
!    integer :: ef2nb_idx1(ne1+1),ef2nb1(2*nf1-nbf1)!,ef2nb_idx2(ne2+1),ef2nb2(2*nf2-nbf2),idx2,ef2nb2mx
!
!
!    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
!    bx = blockDim%x * gridDim%x
!   ! if (e==1) print *,'before', bx
!    !if (e==1) print *,'before', 2*nf2-nbf2,size(ef2nb2),th_id,nsubd,ncs
!
!!if (e ==1) then
!do idxout=e,ncs,bx
!      call get_idx_gpu(index1(idxout),0,g1(idxout),lf1)
!      idx=ef2nb_idx1(g1(idxout))+lf1-1
!
!!      if (idxout==ncs) then
!!        print *, '1', idx,g1,lf1
!!      end if
!
!      call get_idx_gpu(ef2nb1(idx),0,gnb1(idxout),tmp)
!
!
!end do
!end subroutine

!attributes(global) subroutine update_halos_kernel_side2(ef2nb_idx1,ef2nb1, phic1,ef2nb_idx2,ef2nb2, phic2,index1,index2,ncs,ne1,nf1,nbf1,ne2,nf2,nbf2,sz_id1,sz_id2)
!
!    integer :: idxout,g1,g2,lf1,lf2,gnb1,gnb2,tmp,idx
!    integer :: e,bx
!    integer, value :: ne1,nf1,nbf1,ne2,nf2,nbf2,ncs,sz_id1,sz_id2
!    integer :: index1(sz_id1),index2(sz_id2)
!    real :: phic1(ne1+nbf1), phic2(ne2+nbf2)
!    integer :: ef2nb_idx1(ne1+1),ef2nb1(2*nf1-nbf1),ef2nb_idx2(ne2+1),ef2nb2(2*nf2-nbf2),idx2,ef2nb2mx
!
!
!    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
!    bx = blockDim%x * gridDim%x
!   ! if (e==1) print *,'before', bx
!    !if (e==1) print *,'before', 2*nf2-nbf2,size(ef2nb2),th_id,nsubd,ncs
!
!!if (e ==1) then
!do idxout=e,ncs,bx
!      call get_idx_gpu(index1(idxout),0,g1,lf1)
!      idx=ef2nb_idx1(g1)+lf1-1
!
!!      if (idxout==ncs) then
!!        print *, '1', idx,g1,lf1
!!      end if
!
!      call get_idx_gpu(ef2nb1(idx),0,gnb1,tmp)
!
!!      if (idxout==ncs) then
!!        print *, '2', gnb1,tmp
!!      end if
!
!      call get_idx_gpu(index2(idxout),0,g2,lf2)
!      idx=ef2nb_idx2(g2)+lf2-1
!
!!      if (idxout==ncs) then
!!        print *, '3', index2(idxout),idx,g2,lf2
!!      end if
!
!      call get_idx_gpu(ef2nb2(idx),0,gnb2,tmp)
!
!!      if (idxout==ncs) then
!!        print *, '4', gnb2,tmp
!!      end if
!
!      phic1(gnb1)=phic2(g2)
!      phic2(gnb2)=phic1(g1)
!
!
!end do
!    end subroutine

    attributes(global) subroutine update_halos_kernel(ef2nb_idx1,ef2nb1, phic1,ef2nb_idx2,ef2nb2, phic2,index1,index2,ncs,ne1,nf1,nbf1,ne2,nf2,nbf2,sz_id1,sz_id2)

    integer :: idxout,g1,g2,lf1,lf2,gnb1,gnb2,tmp,idx
    integer :: e,bx
    integer, value :: ne1,nf1,nbf1,ne2,nf2,nbf2,ncs,sz_id1,sz_id2
    integer :: index1(sz_id1),index2(sz_id2)
    real :: phic1(ne1+nbf1), phic2(ne2+nbf2)
    integer :: ef2nb_idx1(ne1+1),ef2nb1(2*nf1-nbf1),ef2nb_idx2(ne2+1),ef2nb2(2*nf2-nbf2),idx2,ef2nb2mx


    e = threadIdx%x+(blockIdx%x-1)*blockDim%x
    bx = blockDim%x * gridDim%x
   ! if (e==1) print *,'before', bx
    !if (e==1) print *,'before', 2*nf2-nbf2,size(ef2nb2),th_id,nsubd,ncs

!if (e ==1) then
do idxout=e,ncs,bx
      call get_idx_gpu(index1(idxout),0,g1,lf1)
      idx=ef2nb_idx1(g1)+lf1-1

!      if (idxout==ncs) then
!        print *, '1', idx,g1,lf1
!      end if

      call get_idx_gpu(ef2nb1(idx),0,gnb1,tmp)

!      if (idxout==ncs) then
!        print *, '2', gnb1,tmp
!      end if

      call get_idx_gpu(index2(idxout),0,g2,lf2)
      idx=ef2nb_idx2(g2)+lf2-1

!      if (idxout==ncs) then
!        print *, '3', index2(idxout),idx,g2,lf2
!      end if

      call get_idx_gpu(ef2nb2(idx),0,gnb2,tmp)

!      if (idxout==ncs) then
!        print *, '4', gnb2,tmp
!      end if

!      phic1(gnb1)=phic2(g2)
!      phic2(gnb2)=phic1(g1)


end do
    end subroutine

end module
