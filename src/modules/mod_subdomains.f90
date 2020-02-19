module mod_subdomains
  use mod_cell
  use mod_util
!$ifdef CUDA
  use cudafor
!$endif
  implicit none
  integer, parameter :: CPU=1,GPU=2

  type :: subdomain_t
    integer :: id
    integer :: ne,nf,nbf
    integer, allocatable :: ef2nb(:),ef2nb_idx(:)
    real, allocatable, dimension(:) :: ap,anb,b,phic
    integer :: ARCH=CPU
!$ifdef CUDA
    real, device, allocatable, dimension(:) :: ap_d,anb_d,b_d,phic_d,phic0_d,res_d
    integer, device, allocatable :: ef2nb_d(:),ef2nb_idx_d(:)
    integer :: device_num,stream=0
!$endif
  end type
  type :: intf_t
    integer :: c1,c2
    integer, pointer :: index1(:),index2(:)
    integer :: ncs=0
  end type
 contains
  subroutine construct_subdomains(subdomain,n_subdomains,intf,geom)
    use mod_util
    implicit none
    type(geometry_t) :: geom
    integer :: n_subdomains
    type(intf_t), allocatable :: intf(:,:),intf_tmp(:,:)
    type(subdomain_t), allocatable :: subdomain(:)
    integer :: ef2nb_tmp(2*geom%nf-geom%nbf)
    integer :: g,gf,c,nintf,i,ne,gfnb,lf,l,idx,n,nf,nbf,ncs,idx2,nl,f,cnb,lfnb,cs

    allocate(subdomain(n_subdomains))
    allocate(intf(n_subdomains,n_subdomains))! diagonal is c2b, off diag are intf
    allocate(intf_tmp(n_subdomains,n_subdomains))

    associate(g2gf=>geom%mg%g2gf(1)%p,g2gf_idx=>geom%mg%g2gf(1)%idx,gf2g=>geom%mg%gf2g(1)%p)
!$omp do
    do c=1,n_subdomains
      subdomain(c)%ne=g2gf_idx(c+1)-g2gf_idx(c)
      ne=subdomain(c)%ne
      print *,'ncv in subdomain ',c,ne
      allocate(subdomain(c)%ap(ne))
      allocate(subdomain(c)%b(ne))

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
      allocate(subdomain(c)%phic(ne+nbf))
      allocate(subdomain(c)%ef2nb(2*nf-nbf))
      allocate(subdomain(c)%ef2nb_idx(ne+1))

      ! allocate interface arrays
      do cnb=1,n_subdomains
        ncs=intf(c,cnb)%ncs
        intf(c,cnb)%c1=c
        intf(c,cnb)%c2=cnb

        allocate(intf(c,cnb)%index1(ncs))
        if(c/=cnb) allocate(intf_tmp(c,cnb)%index1(ncs))
      end do
    enddo
!$omp enddo
!$omp do
    do c=1,n_subdomains
      do cnb=1,n_subdomains
        if(c/=cnb) intf(c,cnb)%index2=>intf(cnb,c)%index1
        intf(c,cnb)%ncs=0! ncs is going to be used as a counter next and original ncs is going to be recovered at the end
      end do
    end do
!$omp enddo
!$omp do
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
!$omp enddo
!$omp do
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
!$omp enddo
    ! align the index1 and index2 using intf_tmp reference
!$omp do
    do c=1,n_subdomains
      do cnb=1,n_subdomains
        if(c==cnb) cycle
        ncs=intf(c,cnb)%ncs
        call qsort_key(intf_tmp(c,cnb)%index1,intf(c,cnb)%index1,1,ncs)
        deallocate(intf_tmp(c,cnb)%index1)
      end do
    end do
!$omp enddo
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
        subdomain(c)%b(e)=b(g)
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

    do cs=1,intf%ncs
      call get_idx(intf%index1(cs),0,g1,lf1)
      idx=subd1%ef2nb_idx(g1)+lf1-1
      call get_idx(subd1%ef2nb(idx),0,gnb1,tmp)
      call get_idx(intf%index2(cs),0,g2,lf2)
      idx=subd2%ef2nb_idx(g2)+lf2-1
      call get_idx(subd2%ef2nb(idx),0,gnb2,tmp)

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
