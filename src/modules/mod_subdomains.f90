module mod_subdomains
  use mod_cell
  use mod_util
!$ifdef CUDA
  use cudafor
!$endif
  implicit none
  integer, parameter :: CPU=1,GPU=2,GPU_Unified =3

  type :: subdomain_t
    integer :: id
    integer :: ne,nf,nbf
    !real, allocatable, dimension(:) :: ap,anb,b,phic,phic0
    !integer, allocatable :: ef2nb(:),ef2nb_idx(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ARCH=GPU
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! real, allocatable, managed :: rz(:)
    ! pointers to rz for ap, etc...
    real,allocatable :: ap(:),anb(:),b(:)
    real,allocatable,device :: ap_d(:),anb_d(:),b_d(:)
    real,allocatable,managed :: phic(:),phic0(:),r(:)
    real,allocatable,managed :: res(:),res_max(:)

  !  integer,managed, pointer:: ef2nb(:),ef2nb_idx(:)
!
!    real, allocatable, managed :: ap(:),anb(:),b(:),phic(:),phic0(:),r(:)
    integer, managed,allocatable :: ef2nb(:),ef2nb_idx(:)
 !   integer, allocatable,device :: ef2nb_d(:),ef2nb_idx_d(:)
   ! integer, allocatable,managed :: ef2nb_d(:),ef2nb_idx_d(:)

    integer :: device_num=0
    integer(kind=cuda_Stream_Kind) :: stream,stream1
    integer :: prefetch=0
    integer :: advice=0


  end type


  type :: intf_t
    integer :: c1,c2
    integer, pointer :: index1(:),index2(:)
    integer, managed,allocatable :: index1_d(:),index2_d(:)

  !  integer, device,allocatable :: index1_d(:),index2_d(:),g1_d(:),gnb1_d(:),g2_d(:),gnb2_d(:)

    integer :: ncs=0
  end type
 contains
  subroutine construct_subdomains(subdomain,n_subdomains,intf,geom)
    use mod_util
    use cudafor
    implicit none
    type(geometry_t) :: geom
    integer :: n_subdomains
    type(intf_t), allocatable :: intf(:,:),intf_tmp(:,:)
    type(subdomain_t), allocatable :: subdomain(:)
    integer :: ef2nb_tmp(2*geom%nf-geom%nbf)
    integer :: g,gf,c,nintf,i,ne,gfnb,lf,l,idx,n,nf,nbf,ncs,idx2,nl,f,cnb,lfnb,cs,istat,psize,length,devid

    allocate(subdomain(n_subdomains))

    allocate(intf(n_subdomains,n_subdomains))! diagonal is c2b, off diag are intf
    allocate(intf_tmp(n_subdomains,n_subdomains))

    associate(g2gf=>geom%mg%g2gf(1)%p,g2gf_idx=>geom%mg%g2gf(1)%idx,gf2g=>geom%mg%gf2g(1)%p)

    do c=1,n_subdomains
      subdomain(c)%ne=g2gf_idx(c+1)-g2gf_idx(c)
      ne=subdomain(c)%ne



    allocate(subdomain(c)%ap(ne))
    allocate(subdomain(c)%b(ne))

    subdomain(c)%device_num=c-1


      ! find out the number of control surfaces in interfaces connected to this subdomain
      nf=0
      nbf=0
      do n=g2gf_idx(c),g2gf_idx(c+1)-1
        g=n-g2gf_idx(c)+1
        gf=g2gf(n)
        do idx=geom%ef2nb_idx(gf),geom%ef2nb_idx(gf+1)-1
          call get_idx(geom%ef2nb(idx,1),0,gfnb,lfnb)
          cnb=0
          if(lfnb/=0) cnb=gf2g(gfnb)
          if(cnb==0) then ! this element belongs to a c2b
            intf(c,c)%ncs=intf(c,c)%ncs+1
            nbf=nbf+1
          elseif(c/=cnb) then ! this element belongs to a c2c
            intf(c,cnb)%ncs=intf(c,cnb)%ncs+1
            nbf=nbf+1
          else
          end if
          nf=nf+1
        end do
      end do

      nf=(nf+nbf)/2
      subdomain(c)%nbf=nbf
      subdomain(c)%nf=nf

      allocate(subdomain(c)%anb(2*nf-nbf))
      !allocate(subdomain(c)%phic(ne+nbf))
      !allocate(subdomain(c)%phic0(ne+nbf)) !FOR GPU DEBUG
      allocate(subdomain(c)%ef2nb(2*nf-nbf))
      allocate(subdomain(c)%ef2nb_idx(ne+1))

    !Allocate Unified Memory – accessible from CPU or GPU
    !psize =0
!    istat = cudaMallocManaged(subdomain(c)%rz, 3*ne+2+2*nf-nbf+2*ne+2*nbf+ne, cudaMemAttachGlobal)
!    subdomain(c)%ap => subdomain(c)%rz(1:ne)
!    psize = ne
!    subdomain(c)%b => subdomain(c)%rz(psize+1:psize+ne)
!    psize = psize+ne
!    subdomain(c)%r => subdomain(c)%rz(psize+1:psize+ne)
!    psize = psize+ne
!!    subdomain(c)%res=> subdomain(c)%rz(psize+1:psize+1)
!!    psize = psize+1
!!    subdomain(c)%res_max=> subdomain(c)%rz(psize+1:psize+1)
!!    psize = psize+1
!    subdomain(c)%anb=> subdomain(c)%rz(psize+1:psize+2*nf-nbf)
!    psize = psize+2*nf-nbf
!    subdomain(c)%phic=> subdomain(c)%rz(psize+1:psize+ne+nbf)
!    psize = psize+ne+nbf
!    subdomain(c)%phic0=> subdomain(c)%rz(psize+1:psize+ne+nbf)
!    psize = psize+ne+nbf
!!    subdomain(c)%ef2nb=> subdomain(c)%rz(psize+1:psize+2*nf-nbf)
!!    psize = psize+2*nf-nbf
!!    subdomain(c)%ef2nb_idx=> subdomain(c)%rz(psize+1:psize+ne+1)

!    istat = cudaMallocManaged(subdomain(c)%ap, ne, cudaMemAttachGlobal)
!    istat = cudaMallocManaged(subdomain(c)%b, ne, cudaMemAttachGlobal)
    istat = cudaMallocManaged(subdomain(c)%r, ne, cudaMemAttachHost)
    istat = cudaMallocManaged(subdomain(c)%res, 1, cudaMemAttachHost)
    istat = cudaMallocManaged(subdomain(c)%res_max, 1, cudaMemAttachHost)
!    istat = cudaMallocManaged(subdomain(c)%anb, 2*nf-nbf, cudaMemAttachGlobal)
    istat = cudaMallocManaged(subdomain(c)%phic, ne+nbf, cudaMemAttachHost)
    istat = cudaMallocManaged(subdomain(c)%phic0, ne+nbf, cudaMemAttachHost)
    istat = cudaMallocManaged(subdomain(c)%ef2nb, 2*nf-nbf, cudaMemAttachGlobal)
    istat = cudaMallocManaged(subdomain(c)%ef2nb_idx, ne+1, cudaMemAttachGlobal)
    print *, 'Unified Memory On'
    !istat = cudaMemAdvise(subdomain(c)%rc, 4*ne+3+4*nf-2*nbf+2*ne+2*nbf+ne, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
    do devid=1,4
               ! devid=1

    istat = cudaMemAdvise(subdomain(c)%ef2nb, 2*nf-nbf, cudaMemAdviseSetReadMostly,devid-1)
    istat = cudaMemAdvise(subdomain(c)%ef2nb_idx, ne+1, cudaMemAdviseSetReadMostly,devid-1)
    end do

    istat = cudaMemAdvise(subdomain(c)%r, ne, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)

    if (subdomain(c)%advice==1) then


!    istat = cudaMemAdvise(subdomain(c)%ap, ne, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%b, ne, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%anb, 2*nf-nbf,cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%ef2nb, 2*nf-nbf, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%ef2nb_idx, ne+1, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
    istat = cudaMemAdvise(subdomain(c)%r, ne, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)
!    istat = cudaMemAdvise(subdomain(c)%res, 1, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)
 !   istat = cudaMemAdvise(subdomain(c)%res_max, 1, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)
    istat = cudaMemAdvise(subdomain(c)%phic, ne+nbf, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)
    istat = cudaMemAdvise(subdomain(c)%phic0, ne+nbf, cudaMemAdviseSetPreferredLocation,subdomain(c)%device_num)
!    istat = cudaMemAdvise(subdomain(c)%phic, ne+nbf, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%phic0, ne+nbf, cudaMemAdviseSetPreferredLocation,cudaCpuDeviceId)

!    istat = cudaMemAdvise(subdomain(c)%ap, ne, cudaMemAdviseSetReadMostly,subdomain(c)%device_num)
!    istat = cudaMemAdvise(subdomain(c)%b, ne, cudaMemAdviseSetReadMostly,subdomain(c)%device_num)
!    istat = cudaMemAdvise(subdomain(c)%anb, 2*nf-nbf,cudaMemAdviseSetReadMostly,subdomain(c)%device_num)

   ! istat = cudaMemAdvise(subdomain(c)%r, ne, cudaMemAdviseSetReadMostly,cudaCpuDeviceId)
!    istat = cudaMemAdvise(subdomain(c)%res, 1, cudaMemAdviseSetReadMostly,cudaCpuDeviceId)
 !   istat = cudaMemAdvise(subdomain(c)%res_max, 1, cudaMemAdviseSetReadMostly,cudaCpuDeviceId)
    istat = cudaMemAdvise(subdomain(c)%phic, ne+nbf, cudaMemAdviseSetAccessedBy,cudaCpuDeviceId)
    istat = cudaMemAdvise(subdomain(c)%phic0, ne+nbf, cudaMemAdviseSetAccessedBy,cudaCpuDeviceId)

!    istat = cudaMemAdvise(subdomain(c)%phic, ne+nbf, cudaMemAdviseSetAccessedBy,subdomain(c)%device_num)
!    istat = cudaMemAdvise(subdomain(c)%phic0, ne+nbf, cudaMemAdviseSetAccessedBy,subdomain(c)%device_num)
    end if




      ! allocate interface arrays
      do cnb=1,n_subdomains
        ncs=intf(c,cnb)%ncs
        intf(c,cnb)%c1=c
        intf(c,cnb)%c2=cnb
        allocate(intf(c,cnb)%index1(ncs))
        allocate(intf(c,cnb)%index2(ncs))!!!!! NEW !!!
        length=size(intf(c,cnb)%index1)
      !  print *,'hjk', c,cnb, length
        if(c/=cnb) allocate(intf_tmp(c,cnb)%index1(ncs))
      end do
    enddo

    do c=1,n_subdomains
      do cnb=1,n_subdomains
      !  if(c/=cnb) intf(c,cnb)%index2=intf(cnb,c)%index1
        if(c/=cnb) intf(c,cnb)%index2=>intf(cnb,c)%index1
        intf(c,cnb)%ncs=0! ncs is going to be used as a counter next and original ncs is going to be recovered at the end
      end do
    end do
    !REMOVE POINTER

    do c=1,n_subdomains
      ! fill the connectivity info now
      nbf=0
      ne=subdomain(c)%ne
      subdomain(c)%ef2nb_idx(1)=1
      do n=g2gf_idx(c),g2gf_idx(c+1)-1
        g=n-g2gf_idx(c)+1
        gf=g2gf(n)
        nl=geom%ef2nb_idx(gf+1)-geom%ef2nb_idx(gf)
        subdomain(c)%ef2nb_idx(g+1)=subdomain(c)%ef2nb_idx(g)+nl

        do idx=geom%ef2nb_idx(gf),geom%ef2nb_idx(gf+1)-1
          lf=idx-geom%ef2nb_idx(gf)+1

          call get_idx(geom%ef2nb(idx,1),0,gfnb,lfnb)
          cnb=0
          if(lfnb/=0) cnb=gf2g(gfnb)
          if(c==cnb) then
            idx2=geom%ef2nb_idx(gfnb)+lfnb-1
            ef2nb_tmp(idx2)=index_t(g,lf,0)
          else
            nbf=nbf+1
            ef2nb_tmp(idx)=index_t(ne+nbf,0,0)
          end if
        enddo
      enddo
    enddo

    do c=1,n_subdomains
      do n=g2gf_idx(c),g2gf_idx(c+1)-1
        g=n-g2gf_idx(c)+1
        gf=g2gf(n)

        do idx=geom%ef2nb_idx(gf),geom%ef2nb_idx(gf+1)-1
          lf=idx-geom%ef2nb_idx(gf)+1

          call get_idx(geom%ef2nb(idx,1),0,gfnb,lfnb)
          cnb=0
          if(lfnb/=0) cnb=gf2g(gfnb)

          idx2=subdomain(c)%ef2nb_idx(g)+lf-1
          subdomain(c)%ef2nb(idx2)=ef2nb_tmp(idx)

          if(cnb==0) then ! this element belongs to a c2b
            intf(c,c)%ncs=intf(c,c)%ncs+1
            intf(c,c)%index1(intf(c,c)%ncs)=index_t(g,lf,0)
          elseif(c/=cnb) then ! this element belongs to a c2c
            intf(c,cnb)%ncs=intf(c,cnb)%ncs+1
            cs=intf(c,cnb)%ncs
            intf(c,cnb)%index1(cs)=index_t(g,lf,0)
            if(c<cnb) then
              intf_tmp(c,cnb)%index1(cs)=index_t(gfnb,lfnb,0)
            else
              intf_tmp(c,cnb)%index1(cs)=index_t(gf,lf,0)
            end if
          else ! neighbour is in the same subdomain
          end if
        end do
      end do
    end do
    ! align the index1 and index2 using intf_tmp reference
    do c=1,n_subdomains
      do cnb=1,n_subdomains
        if(c==cnb) cycle
        ncs=intf(c,cnb)%ncs
        call qsort_key(intf_tmp(c,cnb)%index1,intf(c,cnb)%index1,1,ncs)
        deallocate(intf_tmp(c,cnb)%index1)
      end do
    end do

  !   Copy to Device
  !do devid=1,n_subdomains
   !         istat = cudaSetDevice(devid-1)

    do c=1,n_subdomains
      do cnb=c+1,n_subdomains
         !if(c==cnb) cycle
        if (intf(c,cnb)%ncs==0) cycle

        length=size(intf(c,cnb)%index1)
     !   print *,'before', c,cnb, length
!        print *,intf(c,cnb)%index1
       ! allocate(intf(c,cnb)%index1_d(length))
        length=size(intf(c,cnb)%index2)
!                print *,intf(c,cnb)%index2
      !  print *,'after',c,cnb, length
        istat = cudaMallocManaged(intf(c,cnb)%index1_d, length, cudaMemAttachGlobal)
        istat = cudaMallocManaged(intf(c,cnb)%index2_d, length, cudaMemAttachGlobal)

        do devid=1,4

        istat = cudaMemAdvise(intf(c,cnb)%index1_d, length, cudaMemAdviseSetReadMostly,devid-1)
        istat = cudaMemAdvise(intf(c,cnb)%index1_d, length, cudaMemAdviseSetReadMostly,devid-1)
        end do

        intf(c,cnb)%index1_d=intf(c,cnb)%index1
        intf(c,cnb)%index2_d=intf(c,cnb)%index2


     !   allocate(intf(c,cnb)%index2_d(length))

!        istat = cudaMemcpyAsync(intf(c,cnb)%index1_d,intf(c,cnb)%index1,size(intf(c,cnb)%index1),subdomain(c)%stream)
!        istat = cudaMemcpyAsync(intf(c,cnb)%index2_d,intf(c,cnb)%index2,size(intf(c,cnb)%index2),subdomain(c)%stream)
        istat = cudaDeviceSynchronize()
        !print *, istat, c,cnb
      end do
    end do
!end do


    deallocate(intf_tmp)

    end associate
  end subroutine

   subroutine assemble_coef(subdomain,geom,ap,anb,b,phi,nsubd)
    integer :: nsubd
    type(subdomain_t) :: subdomain(nsubd)
    type(geometry_t) :: geom
    real, dimension(*) :: ap,anb,b,phi
    integer :: e,c,idx,idx2,lf,gnb,lfnb,enb,g

    do c=1,nsubd
      do e=1,subdomain(c)%ne
        idx=geom%mg%g2gf(1)%idx(c)+e-1
        g=geom%mg%g2gf(1)%p(idx)
        subdomain(c)%ap(e)=ap(g)
        subdomain(c)%b(e)=b(g) !e = c
        !print *, c,e,g,b(g)
        !pause
        subdomain(c)%phic(e)=phi(g)
        do idx=subdomain(c)%ef2nb_idx(e),subdomain(c)%ef2nb_idx(e+1)-1
          call get_idx(subdomain(c)%ef2nb(idx),0,enb,lfnb)
          lf=idx-subdomain(c)%ef2nb_idx(e)+1
          idx2=geom%ef2nb_idx(g)+lf-1
          subdomain(c)%anb(idx)=anb(idx2)
          if(lfnb==0) then
            call get_idx(geom%ef2nb(idx2,1),0,gnb,lfnb)
            subdomain(c)%phic(enb) = phi(gnb)
          end if
        enddo
      enddo
    enddo

  end subroutine

  subroutine update_halos(intf,subd1,subd2,geom)
    type(subdomain_t) :: subd1,subd2
    type(intf_t) :: intf
    type(geometry_t) :: geom
    integer :: cs,g1,g2,lf1,lf2,gnb1,gnb2,tmp,idx

    if(intf%c1==intf%c2) return
  !  print *,intf%ncs
    do cs=1,intf%ncs
      call get_idx(intf%index1(cs),0,g1,lf1)
      idx=subd1%ef2nb_idx(g1)+lf1-1

!        if (cs==intf%ncs) then
!        print *, '1', idx,g1,lf1
!        end if

      call get_idx(subd1%ef2nb(idx),0,gnb1,tmp)

!              if (cs==intf%ncs) then
!        print *, '2', gnb1,tmp
!        end if

      call get_idx(intf%index2(cs),0,g2,lf2)

      idx=subd2%ef2nb_idx(g2)+lf2-1
!
!         if (cs==intf%ncs) then
!        print *, '3', intf%index2(cs),idx,g2,lf2
!        end if

      call get_idx(subd2%ef2nb(idx),0,gnb2,tmp)

!               if (cs==intf%ncs) then
!        print *, '4', gnb2,tmp
!        end if

      subd1%phic(gnb1)=subd2%phic(g2)
      subd2%phic(gnb2)=subd1%phic(g1)



    end do

    end subroutine

  subroutine update_phi(subdomain,phi,geom,nsubd)
    integer :: nsubd
    type(subdomain_t) :: subdomain(nsubd)
    type(geometry_t) :: geom
    real :: phi(*)
    integer :: c,e,g,idx

    do c=1,nsubd
      do e=1,subdomain(c)%ne
        idx=geom%mg%g2gf(1)%idx(c)+e-1
        g=geom%mg%g2gf(1)%p(idx)
        phi(g) = subdomain(c)%phic(e)
      enddo
    enddo

  end subroutine

end module
