!
!
!  mod_util.f90
!
!  Utility subroutines/functions module
!
module mod_util
    implicit none
    integer, parameter :: sec_lvl=3,num_face_bits=5
!
! Data for unstructured cells
! Element types
    integer :: NODE, BAR_2, BAR_3, TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9
    integer :: TETRA_4, TETRA_10, PYRA_5, PYRA_14,HEXA_8_STR
    integer :: PENTA_6, PENTA_15, PENTA_18, HEXA_8, HEXA_20, HEXA_27
    integer :: MIXED,PYRA_13,NGON_n,NFACE_n
    character (len=32) :: ElementTypeName(0:23)
    parameter (HEXA_8_STR= 1)
    parameter (NODE     =  2)
    parameter (BAR_2    =  3)
    parameter (BAR_3    =  4)
    parameter (TRI_3    =  5)
    parameter (TRI_6    =  6)
    parameter (QUAD_4   =  7)
    parameter (QUAD_8   =  8)
    parameter (QUAD_9   =  9)
    parameter (TETRA_4  = 10)
    parameter (TETRA_10 = 11)
    parameter (PYRA_5   = 12)
    parameter (PYRA_14  = 13)
    parameter (PENTA_6  = 14)
    parameter (PENTA_15 = 15)
    parameter (PENTA_18 = 16)
    parameter (HEXA_8   = 17)
    parameter (HEXA_20  = 18)
    parameter (HEXA_27  = 19)
    parameter (MIXED    = 20)
    parameter (PYRA_13  = 21)
    parameter (NGON_n  = 22)
    parameter (NFACE_n  = 23)
!
    data ElementTypeName / 'Null', 'HEXA_8_STR', &!'UserDefined',                    &
         'NODE', 'BAR_2', 'BAR_3', 'TRI_3', 'TRI_6',                 &
         'QUAD_4', 'QUAD_8', 'QUAD_9', 'TETRA_4', 'TETRA_10',        &
         'PYRA_5', 'PYRA_14', 'PENTA_6', 'PENTA_15', 'PENTA_18',     &
         'HEXA_8', 'HEXA_20', 'HEXA_27', 'MIXED', 'PYRA_13','NGON_n','NFACE_n' /
!
     integer :: ElementTypeDim(0:23)
     data ElementTypeDim /0,3,&
                          0,1,1,2,2,&
                          2,2,2,3,3,&
                          3,3,3,3,3,&
                          3,3,3,3,3,2,3/
! Face vertex data using cgns conventions
     integer :: faces_tetra4(3,4),nfaces_tetra4(4)
     data nfaces_tetra4/3,3,3,3/
     data faces_tetra4/1,3,2, &
                       1,2,4, &
                       2,3,4, &
                       3,1,4 /
!
     integer :: faces_pyra5(4,5),nfaces_pyra5(5)
     data nfaces_pyra5/4,3,3,3,3/
     data faces_pyra5 /1,4,3,2, &
                       1,2,5,0, &
                       2,3,5,0, &
                       3,4,5,0, &
                       4,1,5,0 /
!
     integer :: faces_penta6(4,5),nfaces_penta6(5)
     data nfaces_penta6/4,4,4,3,3/
     data faces_penta6 /1,2,5,4, &
                        2,3,6,5, &
                        3,1,4,6, &
                        1,3,2,0, &
                        4,5,6,0 /
!
     integer :: faces_hexa8(4,6),nfaces_hexa8(6)
     data nfaces_hexa8/4,4,4,4,4,4/
     data faces_hexa8 /1,4,3,2, &
                       1,2,6,5, &
                       2,3,7,6, &
                       3,4,8,7, &
                       1,5,8,4, &
                       5,6,7,8 /
!
     integer :: faces_hexa8_str(4,6) , nfaces_hexa8_str(6),face2halo_str(6)
     data nfaces_hexa8_str/4,4,4,4,4,4/
     data faces_hexa8_str/0,2,3,1, &
                          0,1,5,4, &
                          1,3,7,5, &
                          2,6,7,3, &
                          0,4,6,2, &
                          4,5,7,6/

     data face2halo_str/-4,-2,1,2,-1,4/!b,s,e,n,w,t

     integer :: halo_offset(3,18)
     data halo_offset/ 0, 0,-1,&
                       0,-1, 0,&
                       1, 0, 0,&
                       0, 1, 0,&
                      -1, 0, 0,&
                       0, 0, 1,&
                      -1,-1,0,&
                      -1,+1,0,&
                      1,-1,0,&
                      1,1,0,&
                      -1,0,-1,&
                      -1,0,1,&
                      1,0,-1,&
                      1,0,1,&
                      0,-1,-1,&
                      0,-1,1,&
                      0,1,-1,&
                      0,1,1/
!
     integer :: halo_offset_vx(3,8)
     data halo_offset_vx/ -1, 1,-1,&
                           1, 1,-1,&
                          -1,-1,-1,&
                           1,-1,-1,&
                          -1, 1, 1,&
                           1, 1, 1,&
                          -1,-1, 1,&
                           1,-1, 1/
!
     integer :: halo_offset_i(3,8)
     data halo_offset_i/ 0, 1, 1,&
                         0, 1, 0,&
                         0, 1,-1,&
                         0, 0,-1,&
                         0, 0, 1,&
                         0,-1, 1,&
                         0,-1, 0,&
                         0,-1, -1/
!
     integer :: halo_offset_j(3,8)
     data halo_offset_j/ 1, 0, 1,&
                         0, 0, 1,&
                        -1, 0, 1,&
                         1, 0, 0,&
                        -1, 0, 0,&
                         1, 0,-1,&
                         0, 0,-1,&
                        -1, 0,-1/
!
     integer :: halo_offset_k(3,8)
     data halo_offset_k/ 1, 1, 0,&
                         1, 0, 0,&
                         1,-1, 0,&
                         0,-1, 0,&
                         0, 1, 0,&
                        -1, 1, 0,&
                        -1, 0, 0,&
                        -1,-1, 0/
!
     integer :: faces_tri3(3,1),nfaces_tri3(1)
     data nfaces_tri3/3/
     data faces_tri3 /1,2,3/
!
     integer :: faces_quad4(4,1),nfaces_quad4(1)
     data nfaces_quad4/4/
     data faces_quad4 /1,2,3,4/
!
     integer :: element_nface(20),element_nvx(20)
!                       1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
     data element_nface/6,0,0,0,1,1,1,1,1, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6,0/
     data element_nvx  /8,1,2,3,3,6,4,8,9, 4,10, 5,14, 6,15,18, 8,20,27,0/
!
    type motionsource
      character(len=180) :: modeltype,modelname,modelbc
      integer :: modelindex
      real :: xn(3),yn(3),zn(3)
      real :: xyz(3),uvw(3),uvw0(3),acc(3),ang(3),omega(3),alpha(3),ufar0(3)
      real :: dir(3) ! Direction relative to inertial axis, should always be [1 0 0] for stationary inertial frame
      real :: center(3) ! Center relative to inertial origin, should always be [0 0 0] for stationary inertial frame
      real :: quat(0:3) ! Quaternion for rotation about center relative to inertial axis
      real :: qMatrix(3,3) ! Rotation matrix based on quaternion
    end type motionsource
!
     interface operator(.match.)
        module procedure is_match
     end interface
!
     interface operator(.req.)
        module procedure real_equality
    end interface

    interface operator(.rl.)
        module procedure real_less
    end interface

    interface operator(.in.)
        module procedure is_in_list!, is_in_list_pointer
    end interface

    interface operator(.containsString.)
        module procedure contains_string_fn
    end interface

    interface operator(.veceq.)
        module procedure veceq_int, veceq_real
    end interface

    interface vec_weight
      module procedure vec_weight, vec_weight_proj
    end interface

    interface vec_a_bdotc
      module procedure vec_a_bdotc, vec_a_bdotc_1dim, vec_a_bdotc_1elem
    end interface

    interface int_face1
      module procedure int_face1
    end interface

    interface
      subroutine usleep(useconds) bind(C)
        use iso_c_binding
        implicit none
        integer(c_int32_t), value :: useconds
      end subroutine
    end interface

    type variable
      real, pointer :: ptr(:)
    end type variable

contains
    function real_equality(r1,r2) result(res)
    logical :: res
    real, intent(in) :: r1,r2
    res=.false.
    if(abs(r1-r2) < 1e-6 ) res=.true.
    end function
!
    function real_less(r1,r2) result(res)
    logical :: res
    real, intent(in) :: r1,r2
    res=.false.
    if(r1-r2< 1e-6) res=.true.
    end function
!
    function real_lcm(dt1_i,dt2_i,x1,x2) result(dt) ! real Least Common Multiple
    real :: dt1 , dt2 , dt ,dt1_i,dt2_i
    integer :: x1,x2
    x1=1
    x2=1
    dt1=dt1_i
    dt2=dt2_i
    do while(.not. (dt1 .req. dt2))
      if(dt1 .rl. dt2) then
        x1=x1+1
        dt1=x1*dt1_i
      elseif(dt2 .rl. dt1) then
        x2=x2+1
        dt2=x2*dt2_i
      endif
    enddo
    dt=dt1
    end function

    function contains_string_fn(full_string,candidate_string) result(res)
    implicit none

    character(len=*),intent(in) :: full_string
    character(len=*), intent(in) :: candidate_string
    logical :: res

    res = index(trim(full_string),trim(candidate_string)) > 0

    end function

    function is_in_list_pointer(name,list) result(res)
    logical :: res
    character(len=*) :: name
    character(len=*) , pointer :: list(:)
    integer :: i
    res=.false.
    do i=1,size(list)
      if(trim(name)==trim(list(i)) ) then
        res= .true.
        exit
      endif
    enddo
    end function

    function is_in_list(name,list) result(res)
    logical :: res
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: list(:)
    integer :: i
    res=.false.
    do i=1,size(list)
      if(trim(name)==trim(list(i)) ) then
        res= .true.
        exit
      endif
    enddo
    end function

! This function replaces .match. the operator for real vertices
    function is_match(a,b) result(match)
        real , intent(in) :: a(3),b(3)
        logical :: match
        integer :: i,cnt,n
        real :: r,tolr

        match=.false.
        !tolr=1e-7
        tolr=1e-15

        cnt=0
        do i=1,3
           n = int(log10(abs(b(i)) + tolr))
           !n = max(-7,n-6)
           n = max(-15,n-14)
           tolr = 10.0**n
           if (abs(b(i)-a(i)) < tolr) then
             cnt=cnt+1
           else
             exit
           endif
        enddo
        if(cnt==3) match=.true.
    end function is_match

! used to find if any component of a node matches
    function is_match_point(a,b) result(match)
        real , intent(in) :: a(3),b(3)
        logical :: match(3)
        integer :: i,n
        real :: r,tolr

        match=.false.
        tolr=1e-15

        do i=1,3
           n = int(log10(abs(b(i)) + tolr))
           n = max(-7,n-6)
           tolr = 10.0**n
           if (abs(b(i)-a(i)) < tolr) match(i)=.true.
        enddo

    end function is_match_point

    function getexpnt(r) result(n)
        integer :: n
        real :: r,tolr
          tolr=tiny(r)
          n = int(log10(abs(r) + tolr))
    end function

    function ijk2n(i,j,k,icvs,jcvs) result(ncv)
      integer :: i,j,k,icvs,jcvs,kcvs,ncv

      ncv = i + (j-1)*icvs + (k-1)*icvs*jcvs

    end function

    function ijk2mhalo(i,j,k,icvs,jcvs) result(ncv)
      integer :: i,j,k,icvs,jcvs,kcvs,ncv

      ncv = i+1 + j*(icvs+2) + k*(icvs+2)*(jcvs+2)

    end function

    subroutine n2ijk(n,i,j,k,icvs,jcvs)
        integer :: n,i,j,k,icvs,jcvs

        k=(n-1)/(icvs*jcvs)+1;i=n-(k-1)*icvs*jcvs
        j=(i-1)/icvs+1
        i=mod(i-1,icvs)+1
    end subroutine

    function n2mhalo(n,icvs,jcvs) result(m)
      integer :: n,m,icvs,jcvs
      integer :: i,j,k

      k = (n-1)/(icvs*jcvs)+1; i = n - (k-1)*icvs*jcvs
      j = (i-1)/icvs+1
      i = mod(i-1,icvs)+1

      m = i+1 + j*(icvs+2) + k*(icvs+2)*(jcvs+2)
    end function

    function qijk2n(q,i,j,k,icvs,jcvs) result(ncv)
        integer :: q,i,j,k,icvs,jcvs,kcvs,ncv

        ncv=q+(i-1)*3+(j-1)*3*icvs+(k-1)*3*icvs*jcvs
    end function

    subroutine n2qijk(n,q,i,j,k,icvs,jcvs)
        integer :: n,q,i,j,k,icvs,jcvs

        k=(n-1)/(3*icvs*jcvs)+1;i=n-(k-1)*3*icvs*jcvs
        j=(i-1)/(3*icvs)+1;q=i-(j-1)*3*icvs
        i=(q-1)/3+1
        q=mod(q-1,3)+1
    end subroutine

    function pqijk2n(p,q,i,j,k,icvs,jcvs) result(ncv)
        integer :: p,q,i,j,k,icvs,jcvs,kcvs,ncv

        ncv=p+(q-1)*3+(i-1)*9+(j-1)*9*icvs+(k-1)*9*icvs*jcvs
    end function

    subroutine n2pqijk(n,p,q,i,j,k,icvs,jcvs)
        integer :: n,p,q,i,j,k,icvs,jcvs

        k=(n-1)/(9*icvs*jcvs)+1;i=n-(k-1)*9*icvs*jcvs
        j=(i-1)/(9*icvs)+1;q=i-(j-1)*9*icvs
        i=(q-1)/9+1;p=q-(i-1)*9
        q=(p-1)/3+1
        p=mod(p-1,3)+1
    end subroutine

!----------------------------------------------
! Simple math utilities
!----------------------------------------------
!
!  set_a_0
!
      subroutine set_a_0(a,n)
      implicit none
      integer :: n
      real :: a(n)
      integer :: i
!
      do i=1,n
        a(i) = 0.0
      enddo
!
      end subroutine set_a_0

      subroutine set_a_i4(a,n,i4)
      implicit none
      integer :: n,i4
      integer :: a(n)
      integer :: i
!
      do i=1,n
        a(i) = i4
      enddo
!
      end subroutine set_a_i4
!
!  set_a_ar
!
      subroutine set_a_ar(a,r,n)
      implicit none
      integer :: n
      real :: a(n),r
      integer :: i
!
      do i=1,n
        a(i) = a(i)*r
      enddo
!
      end subroutine set_a_ar
!
!  set_a_b
!
      subroutine set_a_b(a,b,n)
      implicit none
      integer :: n
      real :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = b(i)
      enddo
!
      end subroutine set_a_b
!
      subroutine set_a_b_i4(a,b,n)
      implicit none
      integer :: n
      integer :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = b(i)
      enddo
!
      end subroutine
!
      subroutine set_a_b_c1(a,b,n)
      implicit none
      integer :: n
      character(len=180) :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = b(i)
      enddo
!
      end subroutine
!
!  set_n_m
!
      subroutine set_n_m(n,m,nn)
      implicit none
      integer :: nn
      integer :: n(nn),m(nn)
      integer :: i
!
      do i=1,nn
        n(i) = m(i)
      enddo
!
      end subroutine set_n_m
!
!  set_adivb
!
      subroutine set_adivb(a,b,n)
      implicit none
      integer :: n
      real :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = a(i)/b(i)
      enddo
!
      end subroutine set_adivb
!
!  set_aplusb
!
      subroutine set_aplusb(a,b,n)
      implicit none
      integer :: n
      real :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = a(i)+b(i)
      enddo
!
      end subroutine set_aplusb
!
!  set_a_r
!
      subroutine set_a_r(a,r,n)
      implicit none
      integer :: n
      real :: a(n),r
      integer :: i
!
      do i=1,n
        a(i) = r
      enddo
!
      end subroutine set_a_r

      subroutine set_ab(a,b,n)
      implicit none
      integer :: n
      real :: a(n),b(n)
      integer :: i
!
      do i=1,n
        a(i) = a(i)*b(i)
      enddo
!
      end subroutine set_ab
!
!  set_r_aveb
!
      subroutine set_r_aveb(r,b,n)
      implicit none
      integer :: n
      real :: b(n),r
      real :: rsum
      integer :: i
!
      rsum = 0.0
      do i=1,n
        rsum = rsum + b(i)
      enddo
      r = rsum/n
!
      end subroutine set_r_aveb
!
!  set_a_vertaveb
!
      subroutine set_a_vertaveb(a,b,ni,nj,nk)
      implicit none
      integer :: ni,nj,nk
      real :: a(0:ni+1,0:nj+1,0:nk+1)
      real :: b(ni+1,nj+1,nk+1)
      integer :: i,j,k

      do k=1,nk
        do j=1,nj
          do i=1,ni
            a(i,j,k) = ( b(i,j,k)  +b(i+1,j,k)  +b(i,j+1,k)  +b(i+1,j+1,k) + &
                         b(i,j,k+1)+b(i+1,j,k+1)+b(i,j+1,k+1)+b(i+1,j+1,k+1) ) /8.0
          enddo
        enddo
      enddo
!
      end subroutine set_a_vertaveb
!
!  set_a_swapb
!
      subroutine set_a_swapb(a,b)
      implicit none
      integer :: a,b,s
!
        s = a
        a = b
        b = s
!
      end subroutine set_a_swapb
!
!----------------------------------------------
! Vector related utilities
!----------------------------------------------
      function veceq_int(a,b) result(equals)
        implicit none
        integer, intent(in), dimension(:) :: a,b
        logical :: equals
        integer :: i

        equals = .false.

        if(size(a) /= size(b)) return

        do i = 1,size(a)
          if(a(i) /= b(i)) return
        enddo

        equals = .true.
      end function veceq_int

      function veceq_real(a,b) result(equals)
        implicit none
        real, intent(in), dimension(:) :: a,b
        logical :: equals
        integer :: i

        equals = .false.

        if(size(a) /= size(b)) return

        do i = 1,size(a)
          if(a(i) /= b(i)) return
        enddo

        equals = .true.
      end function
!
!  vec_a_bdotc/vecd_a_bdotc
!
      subroutine vec_a_bdotc(a,b,c,n)
      implicit none
      integer :: n
      real :: a(n),b(3,n),c(3,n)
      integer :: i
!
      do i=1,n
        a(i) = b(1,i)*c(1,i)+b(2,i)*c(2,i) + b(3,i)*c(3,i)
      enddo
!
      end subroutine vec_a_bdotc
!
      !Vector dot product when size(a) = 1
      subroutine vec_a_bdotc_1dim(a,b,c,n)
      implicit none
      real :: a(1),b(3),c(3)
      integer :: n
      !real :: a(1),b(3),c(3)

      a(1) = b(1)*c(1) + b(2)*c(2) + b(3)*c(3)

      end subroutine

      subroutine vec_a_bdotc_1elem(a,b,c,n)
      implicit none
      real :: a, b(3), c(3)
      integer :: n
      !real :: a(1),b(3),c(3)

      a = b(1)*c(1) + b(2)*c(2) + b(3)*c(3)

      end subroutine
!
!  vec_a_bcrossc/vecd_a_bcrossc
!
      subroutine vec_a_bcrossc(a,b,c,n)
      implicit none
      integer :: n
      real :: a(3,n),b(3,n),c(3,n)
      integer :: i
!
      do i=1,n
        a(1,i) = b(2,i)*c(3,i) - b(3,i)*c(2,i)
        a(2,i) = b(3,i)*c(1,i) - b(1,i)*c(3,i)
        a(3,i) = b(1,i)*c(2,i) - b(2,i)*c(1,i)
      enddo
!
      end subroutine vec_a_bcrossc
!  vec_a_magb/vecd_a_magb
!
      subroutine vec_a_magb(a,b,n)
      implicit none
      integer :: n
      real :: a(n),b(3,n)
      integer :: i
!
      do i=1,n
        a(i) = sqrt( b(1,i)*b(1,i)+b(2,i)*b(2,i)+b(3,i)*b(3,i) )
      enddo

      end subroutine vec_a_magb
!
! Function: vecmag
!
      real function vecmag(b)
      implicit none
      real :: b(3)
      vecmag = sqrt( b(1)*b(1)+b(2)*b(2)+b(3)*b(3) )
      end function vecmag
!
! vec_weight
!
      subroutine vec_weight(wt,r0,r1,r2)
      implicit none
      real :: wt,r0(3),r1(3),r2(3)
      real :: ra(3),rb(3),la,lb
!
      ra = r1-r0; rb = r2 - r0
      la = sqrt(ra(1)**2+ra(2)**2+ra(3)**2)
      lb = sqrt(rb(1)**2+rb(2)**2+rb(3)**2)
      if(la+lb>0.) then
        wt = la/(la+lb)
      else
        wt= 0.5
      endif
!
      end subroutine vec_weight
!
      subroutine vec_weight_proj(wt,r0,r1,r2,norm)
      implicit none
      real :: wt,r0(3),r1(3),r2(3),norm(3)
      real :: ra(3),rb(3),la,lb
!
      ra = r1-r0; rb = r2 - r0
      la = abs(dot_product(ra,norm))
      lb = abs(dot_product(rb,norm))
      if(la+lb>0.) then
        wt = la/(la+lb)
      else
        wt= 0.5
      endif
!
      end subroutine vec_weight_proj
!
      subroutine vec_weight_proj3(wt,r0,r1,r2,xn,yn,zn)
      implicit none
      real :: wt(3),r0(3),r1(3),r2(3),xn(3),yn(3),zn(3)
      real :: ra(3),rb(3),la,lb
!
      ra = r1-r0; rb = r2 - r0
      la = abs(dot_product(ra,xn))
      lb = abs(dot_product(rb,xn))
      if(la+lb>0.) then
        wt(1) = la/(la+lb)
      else
        wt(1)= 0.5
      endif
!
      la = abs(dot_product(ra,yn))
      lb = abs(dot_product(rb,yn))
      if(la+lb>0.) then
        wt(2) = la/(la+lb)
      else
        wt(2)= 0.5
      endif
!
      la = abs(dot_product(ra,zn))
      lb = abs(dot_product(rb,zn))
      if(la+lb>0.) then
        wt(3) = la/(la+lb)
      else
        wt(3)= 0.5
      endif
      end subroutine
! gradient at the integration point
      subroutine grad_face(gip,gp,gnb,phip,phinb,rip,rp,rnb)
      implicit none
      real :: gp(3),gnb(3),gip(3),rip(3),rp(3),rnb(3),phip,phinb
      real :: wt,gbar(3),rnorm(3),del

      del = sqrt((rnb(1)-rp(1))**2+(rnb(2)-rp(2))**2+(rnb(3)-rp(3))**2)
      rnorm = (rnb-rp)/del
      call vec_weight(wt,rip,rnb,rp)
      gbar = wt*gp+(1.-wt)*gnb
      gip = gbar+( (phinb-phip)/del-dot_product(gbar,rnorm))*rnorm

      end subroutine
!
      function check_point_in_cyl(xp,yp,zp,p1,p2,radius_sq,length_sq,vec) result(lintersect)
       implicit none
       real :: xp,yp,zp,radius_sq,length_sq,p1(3),p2(3),p(3),vec(3),dist_sq,dot(1),a(3)
       logical :: lintersect

       lintersect = .false.

       p(1) = xp; p(2) = yp; p(3) = zp

       a = p - p1
       call vec_a_bdotc(dot,a,vec,1)

       ! Return if point is outside both caps
       if(dot(1) < 0 .or. dot(1) > length_sq) return

       ! Calculate closest squared distance from axis to point
       dist_sq =  dot_product(a,a) - dot(1)*dot(1)/length_sq

       if (dist_sq > radius_sq) return

       lintersect = .true.

      end function check_point_in_cyl
!
!----------------------------------------------
! Matrix transformation and solution related utilities
!----------------------------------------------
       integer function sgn(x)
       implicit none
       integer :: x
!
       sgn = 0
       if (x>=0) then
         sgn=1
       else if (x<0) then
         sgn=-1
       endif
!
       end function sgn

       real function rsgn(x)
       implicit none
       real :: x
!
       rsgn = 0.
       if (x>=0.) then
         rsgn=1.
       else if (x<0.) then
         rsgn=-1.
       endif
!
       end function rsgn
!
! calc_transform
!
       subroutine calc_transform(transf,tvec)
       implicit none
       integer :: transf(3,3),tvec(3)
       integer :: a,b,c
!
       a = tvec(1)
       b = tvec(2)
       c = tvec(3)
!
       transf(1,1) = sgn(a)*del(a,1)
       transf(2,1) = sgn(a)*del(a,2)
       transf(3,1) = sgn(a)*del(a,3)
!
       transf(1,2) = sgn(b)*del(b,1)
       transf(2,2) = sgn(b)*del(b,2)
       transf(3,2) = sgn(b)*del(b,3)
!
       transf(1,3) = sgn(c)*del(c,1)
       transf(2,3) = sgn(c)*del(c,2)
       transf(3,3) = sgn(c)*del(c,3)
!
       contains
!
!-----
       integer function del(x,y)
       implicit none
       integer :: x,y
!
       del = 0
       if (abs(x)==abs(y)) then
         del=1
       endif
!
       end function del
!
       end subroutine calc_transform
!
! calc_transpose
!
       subroutine calc_transpose(transp,transf)
       implicit none
       integer :: transp(3,3),transf(3,3)
!
       transp(1,1) = transf(1,1)
       transp(2,1) = transf(1,2)
       transp(3,1) = transf(1,3)
!
       transp(1,2) = transf(2,1)
       transp(2,2) = transf(2,2)
       transp(3,2) = transf(2,3)
!
       transp(1,3) = transf(3,1)
       transp(2,3) = transf(3,2)
       transp(3,3) = transf(3,3)
!
       end subroutine calc_transpose
!
      subroutine eigen_symm(w,e,a,n)
! Jacobi method for determination of eigenvalues for an NxN SYMMETRIC matrix a.
! Upper off-diagonal of a is destroyed, lower off-diagonal and main diagonal are untouched.
! w = eigenvalues, e = eigenvectors (columns are vectors for eigenvalues in w, in order).
! n = matrix size.
      implicit none
      integer :: n
      real :: e(n,n),a(n,n),w(n)
!
      real :: sd, so, g, t, c, s, z, theta, thresh, h
      integer :: p,q,r,iter
!Initialize e (eigenvector matrix) to identity matrix
      do p=1,n
        do q=1,n
          if (p==q) then
            e(p,q)=1.0
          else
            e(p,q)=0.0
          endif
        enddo
      enddo
!calculate matrix trace
      sd=0.0
      do p=1,n
        w(p)=a(p,p)
        sd=sd+abs(w(p))
      enddo
      sd=sd**2

      do iter=1,50
! Check for convergence
        so=0.0
        do p=1,n
          do q=p+1,n
            so=so+abs(a(p,q))
          enddo
        enddo
        if (so==0.0) then
          exit
        endif
! Set threshold
        if (iter<4) then
          thresh=0.2*so/n**2
        else
          thresh=0.0
	end if
! Sweep upper diagonal
        do p=1,n
          do q=p+1,n
            g=100*abs(a(p,q))
            if (iter>4 .AND. ((abs(w(p))+g)==abs(w(p))) .AND. ((abs(w(q))+g)==abs(w(q)))) then
              a(p,q)=0
            elseif (abs(a(p,q))>thresh) then
! Establish coefficients
              h=w(q)-w(p)
              if ((abs(h)+g)==abs(h)) then
                t=a(p,q)/h
              else
                theta=0.5*h/a(p,q)
                if (theta<0) then
                  t=-1/((1+theta**2)**0.5-theta)
                else
                  t=1/((1+theta**2)**0.5+theta)
                endif
              c=1/((1+t**2)**0.5)
              s=c*t
              z=t*a(p,q)
              endif
! Apply Transformation
              a(p,q)=0
              w(p)=w(p)-z
              w(q)=w(q)+z
              do r=1,(p-1)
                  t=a(r,p)
                  a(r,p)=c*t-s*a(r,q)
                  a(r,q)=s*t+c*a(r,q)
              enddo
              do r=(p+1),(q-1)
                  t=a(p,r)
                  a(p,r)=c*t-s*a(r,q)
                  a(r,q)=s*t+c*a(r,q)
              enddo
              do r=(q+1),n
                  t=A(p,r)
                  a(p,r)=c*t-s*a(q,r)
                  a(q,r)=s*t+c*a(q,r)
              enddo
! update eigenvectors
               do r=1,n
                 t=e(r,p)
                 e(r,p)=c*t-s*e(r,q)
                 e(r,q)=s*t+c*e(r,q)
               enddo
            endif
          enddo
        enddo
      enddo
!
      end subroutine eigen_symm
! Determinant of a 3x3 matrix
      real function det3(a) result(det)
      implicit none
      real :: a(3,3)
!     Employ Cramers rule
      det = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) - a(1,2)*(a(2,1)*a(3,3)-a(3,1)*a(2,3)) + a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
      end function det3
! 3x3 matrix solution
      subroutine solve3x3(a,x,b)
      implicit none
      real :: a(3,3),x(3),b(3)
      real :: deta,detx,dety,detz
      real :: an(3,3),ao(3,3)
      integer :: i,j

      ao=a
      deta = det3(ao)

      ! the matrix a is already normalized so 1e-6 is actually a small number
      if(abs(deta)<1e-6) then
        x=b/(a(1,1)+a(2,2)+a(3,3))
      else
        an = ao
        do i=1,3
          an(i,1) = b(i)
        enddo
        detx = det3(an)

        an=ao
        do i=1,3
          an(i,2) = b(i)
        enddo
        dety = det3(an)

        an=ao
        do i=1,3
          an(i,3) = b(i)
        enddo
        detz = det3(an)

        x(1) = detx/deta;  x(2) = dety/deta; x(3) = detz/deta
      endif

      end subroutine solve3x3
!
!----------------------------------------------
! Interpolation related utilities
!----------------------------------------------
!
!  int_face1 (calculate ip face coordinates)
!
       subroutine int_face1(rip,x,y,z,i,j,k,ni,nj,nk,iflag)
       implicit none
       integer :: i,j,k,ni,nj,nk,lf,di(4),dj(4),dk(4),v,iflag
       real :: rip(3),dr1(3),dr2(3),area1,area2,areavec(3),geom_cntr1(3),geom_cntr2(3)
       real :: x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk)
!
       select case(iflag)
         case(1)
           lf=5
         case(2)
           lf=3
         case(3)
           lf=2
         case(4)
           lf=4
         case(5)
           lf=1
         case(6)
           lf=6
       end select
       v=faces_hexa8_str(1,lf);di(1)=ibits(v,0,1);dj(1)=ibits(v,1,1);dk(1)=ibits(v,2,1)
       v=faces_hexa8_str(2,lf);di(2)=ibits(v,0,1);dj(2)=ibits(v,1,1);dk(2)=ibits(v,2,1)
       v=faces_hexa8_str(3,lf);di(3)=ibits(v,0,1);dj(3)=ibits(v,1,1);dk(3)=ibits(v,2,1)
       v=faces_hexa8_str(4,lf);di(4)=ibits(v,0,1);dj(4)=ibits(v,1,1);dk(4)=ibits(v,2,1)

       dr1(1) = x(i+di(2),j+dj(2),k+dk(2))-x(i+di(1),j+dj(1),k+dk(1))
       dr1(2) = y(i+di(2),j+dj(2),k+dk(2))-y(i+di(1),j+dj(1),k+dk(1))
       dr1(3) = z(i+di(2),j+dj(2),k+dk(2))-z(i+di(1),j+dj(1),k+dk(1))

       dr2(1) = x(i+di(3),j+dj(3),k+dk(3))-x(i+di(1),j+dj(1),k+dk(1))
       dr2(2) = y(i+di(3),j+dj(3),k+dk(3))-y(i+di(1),j+dj(1),k+dk(1))
       dr2(3) = z(i+di(3),j+dj(3),k+dk(3))-z(i+di(1),j+dj(1),k+dk(1))

       call vec_a_bcrossc(areavec,dr1,dr2,1)
       area1=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
       geom_cntr1 = 0.
       do v=1,3
         geom_cntr1=geom_cntr1+[x(i+di(v),j+dj(v),k+dk(v)),y(i+di(v),j+dj(v),k+dk(v)),z(i+di(v),j+dj(v),k+dk(v))]/3.
       enddo

       dr1(1) = x(i+di(4),j+dj(4),k+dk(4))-x(i+di(1),j+dj(1),k+dk(1))
       dr1(2) = y(i+di(4),j+dj(4),k+dk(4))-y(i+di(1),j+dj(1),k+dk(1))
       dr1(3) = z(i+di(4),j+dj(4),k+dk(4))-z(i+di(1),j+dj(1),k+dk(1))

       call vec_a_bcrossc(areavec,dr2,dr1,1)
       area2=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
       geom_cntr2 = 0.
       do v=1,4
         if(v==2) cycle
         geom_cntr2=geom_cntr2+[x(i+di(v),j+dj(v),k+dk(v)),y(i+di(v),j+dj(v),k+dk(v)),z(i+di(v),j+dj(v),k+dk(v))]/3.
       enddo

       rip = (area1*geom_cntr1+area2*geom_cntr2)/(area1+area2)
!
       end subroutine int_face1
! this is equivalant to int_face1 but with the right face numbering convention
       subroutine int_face2(rip,x,y,z,i,j,k,ni,nj,nk,lf)
       implicit none
       integer :: i,j,k,ni,nj,nk,lf,di(4),dj(4),dk(4),v
       real :: rip(3),dr1(3),dr2(3),area1,area2,areavec(3),geom_cntr1(3),geom_cntr2(3)
       real :: x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk)
!
       v=faces_hexa8_str(1,lf);di(1)=ibits(v,0,1);dj(1)=ibits(v,1,1);dk(1)=ibits(v,2,1)
       v=faces_hexa8_str(2,lf);di(2)=ibits(v,0,1);dj(2)=ibits(v,1,1);dk(2)=ibits(v,2,1)
       v=faces_hexa8_str(3,lf);di(3)=ibits(v,0,1);dj(3)=ibits(v,1,1);dk(3)=ibits(v,2,1)
       v=faces_hexa8_str(4,lf);di(4)=ibits(v,0,1);dj(4)=ibits(v,1,1);dk(4)=ibits(v,2,1)

       dr1(1) = x(i+di(2),j+dj(2),k+dk(2))-x(i+di(1),j+dj(1),k+dk(1))
       dr1(2) = y(i+di(2),j+dj(2),k+dk(2))-y(i+di(1),j+dj(1),k+dk(1))
       dr1(3) = z(i+di(2),j+dj(2),k+dk(2))-z(i+di(1),j+dj(1),k+dk(1))

       dr2(1) = x(i+di(3),j+dj(3),k+dk(3))-x(i+di(1),j+dj(1),k+dk(1))
       dr2(2) = y(i+di(3),j+dj(3),k+dk(3))-y(i+di(1),j+dj(1),k+dk(1))
       dr2(3) = z(i+di(3),j+dj(3),k+dk(3))-z(i+di(1),j+dj(1),k+dk(1))

       call vec_a_bcrossc(areavec,dr1,dr2,1)
       area1=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
       geom_cntr1 = 0.
       do v=1,3
         geom_cntr1=geom_cntr1+[x(i+di(v),j+dj(v),k+dk(v)),y(i+di(v),j+dj(v),k+dk(v)),z(i+di(v),j+dj(v),k+dk(v))]/3.
       enddo

       dr1(1) = x(i+di(4),j+dj(4),k+dk(4))-x(i+di(1),j+dj(1),k+dk(1))
       dr1(2) = y(i+di(4),j+dj(4),k+dk(4))-y(i+di(1),j+dj(1),k+dk(1))
       dr1(3) = z(i+di(4),j+dj(4),k+dk(4))-z(i+di(1),j+dj(1),k+dk(1))

       call vec_a_bcrossc(areavec,dr2,dr1,1)
       area2=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
       geom_cntr2 = 0.
       do v=1,4
         if(v==2) cycle
         geom_cntr2=geom_cntr2+[x(i+di(v),j+dj(v),k+dk(v)),y(i+di(v),j+dj(v),k+dk(v)),z(i+di(v),j+dj(v),k+dk(v))]/3.
       enddo

       rip = (area1*geom_cntr1+area2*geom_cntr2)/(area1+area2)
!
       end subroutine int_face2
!
! set_i_c (fill interface data from cell data)
!
       subroutine set_i_c_uns(phi,phint,index1,ncs1)
        implicit none
        integer :: n , ncs1,e,lf,index1(*)
        real :: phint(*),phi(*)

        do n=1,ncs1
	      call get_idx(abs(index1(n)),0,e,lf)
          phint(n)=phi(e)
        enddo
       end subroutine

       subroutine setc_i_c_uns(i,xyzip,rcint,index1,ncs,ef2nb_idx,ef2nb,ne,nf,nbf)
        implicit none
        integer :: ne,nbf,nf,ncs
        integer :: ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)
        integer :: index1(ncs),i
        real :: xyzip(3,nf),rcint(ncs)
        integer :: e,lf,f,n,idx
!i=1,2,3 for x,y,z respectively
        do n=1,ncs
            call get_idx(iabs(index1(n)),0,e,lf)
            idx=ef2nb_idx(e)+lf-1
            f=iabs(ef2nb(idx,2))
            rcint(n)=xyzip(i,f)
        enddo
       end subroutine
!
       subroutine set_halo_coordinates_uns(rp1,rp2,index1,index12,index2,ei1,nelem1,ei2,nelem2,ncs)
         integer :: ei1,ei2,nelem1,nelem2,ncs
         integer :: index1(ncs),index12(ncs),index2(ncs)
         real :: rp1(nelem1,3),rp2(nelem2,3)
         integer :: n1,n2,e1,e2,tmp

         do n1=1,ncs
            call get_idx(iabs(index1(n1)),0,e1,tmp)
            n2=index12(n1)
            call get_idx(iabs(index2(n2)),0,e2,tmp)
            rp1(ei1+n1-1,:)=rp2(e2,:)
            rp2(ei2+n2-1,:)=rp1(e1,:)
         enddo
       end subroutine
!
       subroutine set_halo_coordinates_bndry_uns(rp1,xyzip,index1,ef2nb_idx,ef2nb,ei,ne,nf,nbf,ncs)
         integer :: ei,ne,nf,nbf,ncs
         real :: rp1(ne+nbf,3),xyzip(3,nf)
         integer :: index1(ncs),ef2nb_idx(ne+1),ef2nb(2*nf-nbf,2)

         integer :: n,e,fl,fg,idx

         do n=1,ncs
           call get_idx(iabs(index1(n)),0,e,fl)
           idx=ef2nb_idx(e)+fl-1
           fg=ef2nb(idx,2)
           rp1(ei+n-1,:)=xyzip(:,fg)
         enddo
       end subroutine
!
! set_i_haloc (fill interface data from halo cell data)
!
       subroutine set_i_haloc(cval,phi,phint,list,dimen,icvs,jcvs,kcvs,list_sz)
       implicit none

       integer :: icvs,jcvs,kcvs,list_sz
       integer :: list(list_sz),dimen
       real :: phi(dimen,0:icvs+1,0:jcvs+1,0:kcvs+1)
       real :: phint(dimen,list_sz)
       character (len=*) :: cval
       integer :: i,j,k,e,l,n,v


       if (cval(1:6)=='center') then
         ! list is actually index1
         do n = 1,list_sz
           call get_idx(abs(list(n)),0,e,l)
           call n2ijk(e,i,j,k,icvs,jcvs)
           phint(1:dimen,n) = phi(1:dimen,i,j,k)
         enddo

       elseif (cval(1:6)=='vertex') then
         ! list is actually vxlist
         do n=1,list_sz
           v=list(n)
           call n2ijk(v,i,j,k,icvs,jcvs)! this icvs,jcvs,kcvs are actually icvs+1,etc
           phint(1:dimen,n) = phi(1:dimen,i,j,k)
         enddo

       else
         write(*,*) 'Error in set_i_haloc: unknown data location'
         stop
       endif
!
       end subroutine set_i_haloc

       subroutine set_i_haloc_uns(phi,phint,index1,ncs1,ng,nbs)
       implicit none
       integer :: ncs1,index1(ncs1),ng,nbs,n,g,ls
       real :: phi(ng+nbs)
       real :: phint(ncs1)

       do n = 1,ncs1
           call get_idx(iabs(index1(n)),0,g,ls)
           phint(n) = phi(g)
       enddo
!
       end subroutine
!
! set_c2haloc (copy cell data to cell with halos)
!
       subroutine set_c2haloc(hcell,cell,ni,nj,nk)
       implicit none
       integer :: ni,nj,nk
       real :: hcell(0:ni+1,0:nj+1,0:nk+1)
       real , intent(in) :: cell(ni,nj,nk)

       integer :: i,j,k

       do k = 1,nk
         do j = 1,nj
           do i = 1,ni
             hcell(i,j,k) = cell(i,j,k)
           enddo
         enddo
       enddo

       end subroutine set_c2haloc
!
       subroutine multiply_mixed_halo(c,a,b,n,m,l)
         implicit none

         real, dimension(:), intent(in) :: a,b
         real, dimension(:), intent(out) :: c
         integer, intent(in) :: n,m
         integer, intent(in), optional :: l

         integer :: i,j,k,idx,idx_h

           do idx = 1,n
             c(idx) = a(idx)*b(idx)
           enddo

       end subroutine

       subroutine multiply_accumulate_mixed_halo(c,a,b,n,m,l)
         implicit none

         real, dimension(:), intent(in) :: a,b
         real, dimension(:), intent(inout) :: c
         integer, intent(in) :: n,m
         integer, intent(in), optional :: l

         integer :: i,j,k,idx,idx_h

           do idx = 1,n
             c(idx) = c(idx) + a(idx)*b(idx)
           enddo

       end subroutine

       subroutine divide_mixed_halo(c,a,b,n,m,l)
         implicit none

         real, dimension(:), intent(in) :: a,b
         real, dimension(:), intent(out) :: c
         integer, intent(in) :: n,m
         integer, intent(in), optional :: l

         integer :: i,j,k,idx,idx_h

           do idx = 1,n
             c(idx) = a(idx)/b(idx)
           enddo

       end subroutine

       subroutine subtract_mixed_halo(c,a,b,n,m,l)
         implicit none

         real, dimension(:), intent(in) :: a,b
         real, dimension(:), intent(out) :: c
         integer, intent(in) :: n,m
         integer, intent(in), optional :: l

         integer :: i,j,k,idx,idx_h

           do idx = 1,n
             c(idx) = a(idx)-b(idx)
           enddo

       end subroutine

       subroutine multiply_subtract_mixed_halo(c,a,b,n,m,l)
         implicit none

         real, dimension(:), intent(in) :: a,b
         real, dimension(:), intent(inout) :: c
         integer, intent(in) :: n,m
         integer, intent(in), optional :: l

         integer :: i,j,k,idx,idx_h

           do idx = 1,n
             c(idx) = c(idx) - a(idx)*b(idx)
           enddo

       end subroutine
!
! ----------------------------------
!
     subroutine face_vxlist(vxlist,t,f,lst,nl)
! For an element of type t at current face f, return the vertex list
     implicit none
!
     integer :: t,f
     integer :: vxlist(*)
!
     integer ::nl,lst(*)
!
     integer :: l
!
     if (t==TETRA_4) then
       if (f>4)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_tetra4(f)
       do l=1,nl
         lst(l) = vxlist(faces_tetra4(l,f))
       enddo
     elseif (t==PYRA_5) then
       if (f>5)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_pyra5(f)
       do l=1,nl
         lst(l) = vxlist(faces_pyra5(l,f))
       enddo
     elseif (t==PENTA_6) then
       if (f>5)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_penta6(f)
       do l=1,nl
         lst(l) = vxlist(faces_penta6(l,f))
       enddo
     elseif (t==HEXA_8) then
       if (f>6)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_hexa8(f)
       do l=1,nl
         lst(l) = vxlist(faces_hexa8(l,f))
       enddo
     elseif (t==TRI_3) then
       if (f>1)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_tri3(f)
       do l=1,nl
         lst(l) = vxlist(faces_tri3(l,f))
       enddo
     elseif (t==QUAD_4) then
       if (f>1)then
         write(*,*) 'Error in face number in face_vxlist ...'
       endif
       nl = nfaces_quad4(f)
       do l=1,nl
         lst(l) = vxlist(faces_quad4(l,f))
       enddo
     else
       write(*,*) 'Unsupported element type in face_vxlist ...',t
       stop
     endif
!
     end subroutine face_vxlist

    function index_t(group_no,side_no,lvl)

        integer :: index_t,lvl
        integer , intent(in) :: group_no,side_no

        index_t=0
        index_t=ishft(group_no,min(2*lvl+num_face_bits,16))
        index_t=ior(index_t,side_no)

    end function index_t
!
     subroutine get_idx(idx,lvl,group_no,side_no)

        integer , intent(out) :: group_no,side_no
        integer , intent(in) :: idx
        integer :: lvl

        side_no=iand(idx,2**(min(2*lvl+num_face_bits,16))-1)
        group_no=ishft(idx,-(min(2*lvl+num_face_bits,16)))

     end subroutine get_idx
!
     function face_compare(e1,f1,t1,e2vx1,e2,f2,t2,e2vx2) result(res)
!
! Face 1 node list and face 2 node lis
! Check if the two faces share common nodes
     implicit none
!
     integer :: idx1,idx2,e1,e2,f1,f2
     integer :: nl1,nl2,t1,t2,e2vx1(*),e2vx2(*)
     integer :: lst1(4),lst2(4)
     logical :: res
     integer :: cnt,i,j

    res = .false.
    if(e1==e2) return
    call face_vxlist(e2vx1,t1,f1,lst1,nl1)
    call face_vxlist(e2vx2,t2,f2,lst2,nl2)
!
     if (nl1==nl2) then
       cnt = 0
       do i=1,nl1
         do j=1,nl2
           if (lst1(i)==lst2(j)) then
             cnt = cnt + 1;exit
           endif
         enddo
       enddo
       if (cnt == nl1) res = .true.
     endif
!
     end function face_compare
!----------------------------------------------
! String related utilities
!----------------------------------------------
      subroutine convert_char2int(cold,inew,n)
      implicit none
      integer :: n,inew(n),i,j
      character (len=n) :: cold

        inew = 0
        do i = 1,len_trim(cold)
          inew(i) = ichar(cold(i:i))
        enddo

      end subroutine convert_char2int
!
!----------------------------------------------
      subroutine convert_int2char(iold,cnew,n)
      implicit none
      integer :: n,iold(n),i,j
      character (len=n) :: cnew

        cnew = ''
        do i = 1,n
          if (iold(i) /= 0) then
            cnew(i:i) = char(iold(i))
          elseif (iold(i) == 0) then
            cnew(i:i) = ''
          endif
        enddo

      end subroutine convert_int2char
!
!----------------------------------------------
      subroutine check_casedir(dirname)
      implicit none
      character (len=*) :: dirname
      integer :: n

        do n=1,len_trim(dirname)
          if (dirname(n:n)=='/') call system('mkdir '//dirname(1:n-1)//' 2> /dev/null')
        enddo
        call system('mkdir '//trim(dirname)//' 2> /dev/null')

      end subroutine check_casedir
!
!----------------------------------------------
! Randomization Utilities
!----------------------------------------------
!
    subroutine rand_normal(c,n,m,mean,stdev)
! creates a size n 1D array containing randomly distributed numbers with
! prescribed mean and standard deviation
      implicit none
      real :: mean,stdev
      integer :: i,j,n,m
!
      real :: c(n,m)
      real :: temp(2),r,theta
!
      if(stdev <= 0.0d0) THEN
        write(*,*) "Standard Deviation must be +ve"
        stop
      else
        do j=1,m
          do i=1,n
            call RANDOM_NUMBER(temp)
            r=(-2.0*log(temp(1)))**0.5
            theta = 2.0*3.1415*temp(2)
            c(i,j)= mean+stdev*r*sin(theta)
          enddo
        enddo
      endif
!
    end subroutine rand_normal
!
    subroutine rand_nbr_2D(c,n,m)
! creates a size n 2D array containing uniform randomly distributed numbers
! c = output array, of size nxm
      implicit none
      integer :: i,j,n,m
      real :: c(n,m),temp(1)
!
      do j=1,m
        do i=1,n
          call RANDOM_NUMBER(temp)
          c(i,j)= temp(1)
        enddo
      enddo
!
    end subroutine rand_nbr_2D
!
    subroutine rand_nbr_1D(c,n)
! creates a size n 1D array containing uniform randomly distributed numbers
! c = output array, of size n
      implicit none
      integer :: i,n
      real :: c(n),temp(1)
!
      do i=1,n
        call RANDOM_NUMBER(temp)
        c(i) = temp(1)
      enddo
!
    end subroutine rand_nbr_1D
!
    subroutine iswap(a,b)
      integer :: a,b

      a=a+b
      b=a-b
      a=a-b
    end subroutine

    subroutine rswap(a,b)
      real :: a,b

      a=a+b
      b=a-b
      a=a-b
    end subroutine

  ! sort b against key: quick sort by key
  recursive subroutine qsort_key(key,b,i,f)
    integer :: key(*),b(*)
    integer :: i,f,n,p,itmp

    if(i>=f) return
    n=i
    p=f! pivot
    do while(n<p)
      if(key(n)>key(p)) then
        call iswap(key(n),key(p))
        call iswap(b(n),b(p))
        p=p-1
        if(p>n) then
          call iswap(key(n),key(p))
          call iswap(b(n),b(p))
        endif
      else
        n=n+1
      endif
    enddo
!
    call qsort_key(key,b,i,p-1)
    call qsort_key(key,b,p+1,f)
  end subroutine

  recursive subroutine qsort_rkey(key,b,i,f)
    real :: key(*)
    integer :: b(*)
    integer :: i,f,n,p,itmp

    if(i>=f) return
    n=i
    p=f! pivot
    do while(n<p)
      if(key(n)>key(p)) then
        call rswap(key(n),key(p))
        call iswap(b(n),b(p))
        p=p-1
        if(p>n) then
          call rswap(key(n),key(p))
          call iswap(b(n),b(p))
        endif
      else
        n=n+1
      endif
    enddo
!
    call qsort_rkey(key,b,i,p-1)
    call qsort_rkey(key,b,p+1,f)
  end subroutine

  recursive subroutine qsort_rkey3(key,b1,b2,b3,i,f)
    real :: key(*)
    integer :: b1(*),b2(*),b3(*)
    integer :: i,f,n,p,itmp

    if(i>=f) return
    n=i
    p=f! pivot
    do while(n<p)
      if(key(n)>key(p)) then
        call rswap(key(n),key(p))
        call iswap(b1(n),b1(p))
        call iswap(b2(n),b2(p))
        call iswap(b3(n),b3(p))
        p=p-1
        if(p>n) then
          call rswap(key(n),key(p))
          call iswap(b1(n),b1(p))
          call iswap(b2(n),b2(p))
          call iswap(b3(n),b3(p))
        endif
      else
        n=n+1
      endif
    enddo
!
    call qsort_rkey3(key,b1,b2,b3,i,p-1)
    call qsort_rkey3(key,b1,b2,b3,p+1,f)
  end subroutine

  subroutine qsort_key_nRec(key,b,n)
    integer, parameter :: max_levels=1000
    integer :: n,key(0:n-1),b(0:n-1)
    integer :: beg(0:max_levels-1),end(0:max_levels-1)
    integer :: i,L,R,piv,b0

    i=0
    beg(0)=0;end(0)=n
    do while(i>=0)
      L=beg(i);R=end(i)-1
      if(L<R) then
        piv=key(L)
        b0=b(L)
        if(i==max_levels-1) print *,'wtf!'
        do while(L<R)
          do while(key(R)>=piv .and. L<R)
            R=R-1
          enddo
          if(L<R) then
            key(L)=key(R)
            b(L)=b(R)
            L=L+1
          end if
          do while(key(L)<=piv .and. L<R)
            L=L+1
          enddo
          if(L<R) then
            key(R)=key(L)
            b(R)=b(L)
            R=R-1
          end if
        end do
        key(L)=piv
        b(L)= b0
        beg(i+1)=L+1
        end(i+1)=end(i)
        end(i)=L
        i=i+1
        if(end(i)-beg(i)>end(i-1)-beg(i-1)) then
          call iswap(beg(i),beg(i-1))
          call iswap(end(i),end(i-1))
        end if
      else
        i=i-1
      end if
    enddo

  end subroutine

  subroutine qsort_nRec(arr,n)
    ! Implementation by Darel Rex Finley.http://alienryderflex.com/quicksort/
    integer, parameter :: max_levels=1000
    integer :: n,arr(0:n-1)
    integer :: beg(0:max_levels-1),end(0:max_levels-1)
    integer :: i,L,R,piv

    i=0
    beg(0)=0;end(0)=n
    do while(i>=0)
      L=beg(i);R=end(i)-1
      if(L<R) then
        piv=arr(L)

        do while(L<R)
          do while(arr(R)>=piv .and. L<R)
            R=R-1
          enddo
          if(L<R) then
            arr(L)=arr(R)
            L=L+1
          end if
          do while(arr(L)<=piv .and. L<R)
            L=L+1
          enddo
          if(L<R) then
            arr(R)=arr(L)
            R=R-1
          end if
        end do
        arr(L)=piv
        beg(i+1)=L+1
        end(i+1)=end(i)
        end(i)=L
        i=i+1
        if(end(i)-beg(i)>end(i-1)-beg(i-1)) then
          call iswap(beg(i),beg(i-1))
          call iswap(end(i),end(i-1))
        end if
      else
        i=i-1
      end if
    enddo

  end subroutine

  subroutine insertion_sort_key(key,b,i,f)
    integer :: key(*),b(*)
    integer :: i,f,n,m

    do n=i,f
      do m=n-1,i,-1
        if(key(m)>key(m+1)) then
          call iswap(key(m),key(m+1))
          call iswap(b(m),b(m+1))
        else
          exit
        end if
      end do
    end do
  end subroutine

  recursive subroutine qsort_key2(key,b,c,i,f)
    integer :: key(*),b(*),c(*)
    integer :: i,f,n,p,itmp

    if(i>=f) return
    n=i
    p=f! pivot
    do while(n<p)
      if(key(n)>key(p)) then
        call iswap(key(n),key(p))
        call iswap(b(n),b(p))
        call iswap(c(n),c(p))
        p=p-1
        if(p>n) then
          call iswap(key(n),key(p))
          call iswap(b(n),b(p))
          call iswap(c(n),c(p))
        endif
      else
        n=n+1
      endif
    enddo
!
    call qsort_key2(key,b,c,i,p-1)
    call qsort_key2(key,b,c,p+1,f)
  end subroutine

  subroutine qsort_key_2way(key,b,c,i,f)
    integer :: key(*),b(*),c(*)
    integer :: i,f,n,n0,key0

    call qsort_key2(key,b,c,i,f)
    n0=1
    key0=key(i)-1
    do n=i,f
      if(key(n)/=key0) then
        key0=key(n)
        call qsort_key(b,c,n0,n-1)
        n0=n
      endif
    enddo
    ! last one
    call qsort_key(b,c,n0,n-1)
  end subroutine

  recursive subroutine qsort(a,i,f)
    integer ::a(*)
    integer :: i,f,n,p,itmp

    if(i>=f) return
    n=i
    p=f! pivot
    do while(n<p)
      if(a(n)>a(p)) then
        call iswap(a(n),a(p))
        p=p-1
        if(p>n) call iswap(a(n),a(p))
      else
        n=n+1
      endif
    enddo
!
    call qsort(a,i,p-1)
    call qsort(a,p+1,f)
  end subroutine

! binary search of the value x in a sorted array x. Spits out the index of the element i.e. a(res):=x
  recursive function bsearch(a,x,i,f) result(res)
    integer :: a(*),x,i,f,p,res

    res=0
    if(i>f) return
    p=(f+i)/2
    if(x<a(p)) then
      res=bsearch(a,x,i,p-1)
    elseif(x>a(p)) then
      res=bsearch(a,x,p+1,f)
    else
      res=p
    endif
  end function

  function hf_anb(hf,lvl,cvs,ef2nb_idx,ef2nb) result(idx)
    ! this function returns the position in the anb coefficient array
    integer, intent(in) :: hf,lvl,cvs(3)
    integer :: ef2nb_idx(*),ef2nb(*)
    integer :: e,lf,idx

    call get_idx(hf,lvl,e,lf)
    if(ef2nb_idx(1)/=0) then ! uns
      idx=ef2nb_idx(e)+lf -1
    else
      idx=(lf-1)*cvs(1)*cvs(2)*cvs(3)+e
    endif

  end function

  subroutine hf_info(hf,lvl,cvs,idx,e,lf,gf,enb,lfnb,h,hnb,ef2nb_idx,ef2nb)! half-face info
    integer, intent(in) :: hf,lvl,cvs(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer, intent(out) :: e,h,lf,enb,lfnb,gf,hnb,idx
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    enb = 0
    lfnb = 0
    call get_idx(hf,lvl,e,lf)
! non-zero valued of integers are equivalent to true
      h = e
      idx=ef2nb_idx(e)+lf -1
      i = 2*cvs(2)-cvs(3)+idx
      gf = ef2nb(i)
      call get_idx(ef2nb(idx),lvl,enb,lfnb)
      hnb = enb
      if(lfnb==0) enb = enb - cvs(1)
    end subroutine

  function get_enb(hf,lvl,cvs,e,lf,ef2nb_idx,ef2nb) result(enb)! half-face info
! enb is  the input for the index_inv which returns cs
    integer :: hf,lvl,cvs(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer :: e,lf,enb,lfnb
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    enb = 0
    lfnb = 0! face number for halo elements is 0
    call get_idx(hf,lvl,e,lf)


      i = ef2nb_idx(e)+lf -1
      call get_idx(ef2nb(i),lvl,enb,lfnb)
      if(lfnb==0) enb=enb-cvs(1)

  end function
! get helo
  function get_hnb(hf,lvl,dim,cvs,e,lf,ef2nb_idx,ef2nb) result(hnb)! half-face info
    integer :: hf,lvl,cvs(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer :: e,lf,hnb,dim
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    call get_idx(hf,lvl,e,lf)


      i = ef2nb_idx(e)+lf -1
      call get_idx(ef2nb(i),lvl,hnb,k)
     ! hnb=(dim-1)*(cvs(1)+cvs(3))+hnb ! this is for uvw array
      hnb=(hnb-1)*dim+1   ! this is for vector array, bodyforce or grad_vec

  end function

  function get_hnb_pad(hf,lvl,dim,cvs,pad,e,lf,ef2nb_idx,ef2nb) result(hnb)! half-face info
    integer :: hf,lvl,cvs(3),pad(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer :: e,lf,hnb,dim
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    call get_idx(hf,lvl,e,lf)

      i = ef2nb_idx(e)+lf -1
      call get_idx(ef2nb(i),lvl,hnb,k)
     ! hnb=(dim-1)*(cvs(1)+cvs(3))+hnb
      hnb=(hnb-1)*dim+1

  end function

  function get_h_cvc(hf,lvl,dim,cvs,e,lf,ef2nb_idx,ef2nb) result(h)! half-face info
    integer :: hf,lvl,cvs(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer :: e,lf,h,dim
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    call get_idx(hf,lvl,e,lf)


     ! h=(dim-1)*(cvs(1)+cvs(3))+e
      h=(e-1)*dim+1
  end function

  function get_h_pad(hf,lvl,dim,cvs,pad,e,lf,ef2nb_idx,ef2nb) result(h)! half-face info
    integer :: hf,lvl,cvs(3),pad(3)! half-face index, element type, multigrid level, cell size info uns=(ne,nf,nbf), str=(icvs,jcvs,kcvs)
    integer :: ef2nb_idx(*),ef2nb(*)! null for str
    integer :: e,lf,h,dim
! element number, halo element index,local face number, neighbour element number, local face number of neighbor element, global face number, halo index of neighbour element
    integer :: i,j,k

    call get_idx(hf,lvl,e,lf)

    if(ef2nb_idx(1)/=0) then
     ! h=(dim-1)*(cvs(1)+cvs(3))+e
      h=(e-1)*dim+1
    endif

  end function

  function e2h(e,cvs,ef2nb_idx) result(h)
    integer :: e,h,cvs(3),ef2nb_idx(*)

      h=e

  end function

  function get_ne(cvs,ef2nb_idx) result(ne)
    integer :: cvs(3),ef2nb_idx(*),ne

    if(ef2nb_idx(1)/=0) then
      ne = cvs(1)
    endif
  end function

  function get_nh(cvs,ef2nb_idx) result(nh)
    integer :: cvs(3),ef2nb_idx(*),nh

    if(ef2nb_idx(1)/=0) then
      nh = cvs(1)+cvs(3)
    endif
  end function

  function str(cin) result(cout)
   character(len=*),intent(in) :: cin
   character(len=:),allocatable :: cout

   allocate(cout, source=trim(cin))

  end function
 ! functions for integer flag
 subroutine set_iflag(flag,e)
  integer :: flag(*)
  integer :: idx,pos,e

  idx=(e-1)/(kind(flag)*8)+1
  pos=mod(e-1,kind(flag)*8)
  flag(idx)=ibset(flag(idx),pos)
 endsubroutine

 subroutine clr_iflag(flag,e)
  integer :: flag(*)
  integer :: idx,pos,e

  idx=(e-1)/(kind(flag)*8)+1
  pos=mod(e-1,kind(flag)*8)
  flag(idx)=ibclr(flag(idx),pos)
 endsubroutine

 function get_iflag(flag,e) result(lflag)
  integer :: flag(*)
  integer :: idx,pos,e
  logical :: lflag

  idx=(e-1)/(kind(flag)*8)+1
  pos=mod(e-1,kind(flag)*8)
  lflag= ibits(flag(idx),pos,1)==1
 endfunction

 function transform_FoR(v1,mtn1,mtn2) result(v1_2)
   real, dimension(3) :: v1,v1_2
   type(motionsource) :: mtn1,mtn2
   real, dimension(3) :: dtmp

   dtmp=v1(1)*mtn1%xn+v1(2)*mtn1%yn+v1(3)*mtn1%zn
   v1_2=[dot_product(dtmp,mtn2%xn),dot_product(dtmp,mtn2%yn),dot_product(dtmp,mtn2%zn)]

 end function

end module mod_util

module mod_dll
    implicit none

    type :: value_t
      integer :: id =0
    end type

    type :: node_t
        class(value_t), pointer :: val=>null()
        type(node_t), pointer :: prev=>null()
        type(node_t), pointer :: next =>null()
    end type

    type :: iterator_t
        type(node_t), pointer :: node=>null()
    end type

    interface assignment(=)
      module procedure assign_val
    end interface

    interface node_t
      module procedure construct_node_t
    end interface
! (head)--prev--> .....<--next--(iterator)--prev--> ...<--next--(tail)
    type :: dll_t
        type(node_t), pointer :: head => null()
        type(node_t), pointer :: tail => null()
        integer :: length=0
      contains
        procedure :: get
        procedure :: pop
        procedure :: push
        procedure :: push_front
        procedure :: insert
        procedure :: rewind
        procedure :: iterate
        procedure :: destroy
    endtype

    interface value_t
      module procedure construct_value_t
    end interface
contains
    function construct_value_t(id) result(val)
     type(value_t), pointer :: val
     integer :: id

     allocate(val)
     val%id = id

    end function

    subroutine assign_val(val,iterator)
      class(value_t), pointer, intent(out) :: val
      type(iterator_t), intent(in) :: iterator

      val => iterator%node%val
    end subroutine

    function construct_node_t(val) result(node)
        type(node_t), pointer :: node
        class(value_t), pointer :: val

        allocate(node)
        node%val => val

    end function

    function get(this,id) result(val)
      class(dll_t) :: this
      class(value_t), pointer :: val
      integer :: id
      type(iterator_t) :: iterator

      nullify(val)
      if(this%rewind(iterator)) then
      do
        val = iterator
        if(val%id==id) exit
        if(.not. this%iterate(iterator)) exit
      enddo
      endif
    end function

    logical function pop(this,it,val)
      class(dll_t) :: this
      class(value_t), pointer :: val
      type(iterator_t) :: it
      type(node_t), pointer :: tmp

      pop = this%length>0
      if(pop) then
        val = it
        if(associated(it%node%prev) .and. associated(it%node%next)) then
          it%node%prev%next => it%node%next
          it%node%next%prev => it%node%prev
          tmp => it%node%next
          nullify(it%node%val)
          deallocate(it%node)
          it%node => tmp
        elseif(associated(it%node%prev)) then
          nullify(it%node%prev%next)
          this%head => it%node%prev
          nullify(it%node%val)
          nullify(it%node%prev)
          deallocate(it%node)
          it%node => this%head
        elseif(associated(it%node%next)) then
          nullify(it%node%next%prev)
          this%tail => it%node%next
          nullify(it%node%val)
          deallocate(it%node)
          it%node => this%tail
        else
          nullify(it%node%val)
          deallocate(this%head)
          nullify(it%node)
          nullify(this%head)
          nullify(this%tail)
        endif
        this%length = this%length - 1
      endif

    end function

    subroutine push(this,val)! push_back
      class(dll_t) :: this
      class(value_t), pointer :: val
      type(node_t), pointer :: node

      node => node_t(val)
      if(associated(this%tail)) then
        this%tail%prev=>node
        node%next => this%tail
        this%tail => node
      else
        this%head=>node
        this%tail=>node
      end if
      this%length = this%length+1

    end subroutine

    subroutine push_front(this,val)
      class(dll_t) :: this
      class(value_t), pointer :: val
      type(node_t), pointer :: node

      node => node_t(val)
      if(associated(this%head)) then
        this%head%next=>node
        node%prev => this%head
        this%head => node
      else
        this%head=>node
        this%tail=>node
      end if
      this%length = this%length+1

    end subroutine

    subroutine insert(this,val,it)
      class(dll_t) :: this
      class(value_t), pointer :: val
      type(node_t), pointer :: node
      type(iterator_t) :: it

      if(this%length>0) then
        node => node_t(val)
        if(associated(it%node%prev)) then
          it%node%prev%next => node
          node%prev => it%node%prev
          it%node%prev=>node
          node%next => it%node
        else
          this%tail=>node
          it%node%prev=>node
          node%next => it%node
        endif
      else
        call this%push(val)
        node => this%head
      endif

      it%node => node
      this%length = this%length +1

    end subroutine


    logical function rewind(this,iterator)
      class(dll_t) :: this
      type(iterator_t) :: iterator

      rewind = this%length>0
      if(.not. rewind) return
      iterator%node => this%head

    end function

    logical function iterate(this,it)
      class(dll_t) :: this
      type(iterator_t) :: it

      iterate =.false.
      if(this%length==0) return

      if(associated(it%node%prev)) then
        it%node => it%node%prev
        iterate=.true.
      else
        it%node => this%head
      endif

    end function

    subroutine destroy(this)
      class(dll_t) :: this
      class(value_t), pointer :: val
      type(iterator_t) :: iterator

      nullify(val)
      if(this%rewind(iterator)) then
      do while(this%pop(iterator,val))
        deallocate(val)
      enddo
      endif
    end subroutine

end module
