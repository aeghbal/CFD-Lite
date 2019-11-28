
module mod_cgns
  implicit none
  include 'cgnslib_f.h'

    type cgnsdb
      integer :: base=1
      integer :: unit=1
      character(len=255) :: base_name
      character(len=180) :: file

    end type cgnsdb

    integer :: ier
 contains
!----------------------  cgns_db_open
  subroutine cgns_db_open(cg)
    implicit none
    type (cgnsdb) :: cg
    integer :: n,cell_dim,phys_dim,ndsc,idx
    character (len=180) :: path,text,name

    call cg_open_f(cg%file, CG_MODE_MODIFY, cg%unit, ier)
    if (ier .eq. ERROR) call cg_error_exit_f

! Read in database info (there can only be one)
    call cg_base_read_f(cg%unit, n, cg%base_name, cell_dim, phys_dim, ier)

  end subroutine cgns_db_open

  subroutine cgns_element_info(cg,e2vx,esection,etype,ne2vx,nelem,ne2vxmax,nsec)
!
    implicit none
    integer :: nelem,ne2vxmax,nsec
    integer :: e2vx(ne2vxmax,nelem),esection(2,nsec),etype(nsec),ne2vx(nsec),n
    type (cgnsdb) :: cg
!
!
    integer :: s,i,j,k
    integer, allocatable :: elements(:)
    integer :: eldatasize,type,start,end,nbndry,parent_flag,npe
    integer :: nsections,zone
    character (len=32) :: elsname
!
    e2vx =0
    zone=1
!
!--- zone attribute:  GOTO zone node
    call cg_goto_f(cg%unit, cg%base, ier, 'Zone_t', zone, 'end')
    if (ier .eq. ERROR) call cg_error_exit_f
!
!--- zone element attribute:  GOTO Elements_t node
    call cg_goto_f(cg%unit, cg%base, ier, 'Zone_t', zone, 'Elements_t', 1, 'end')
    if (ier .eq. ERROR) call cg_error_exit_f
!
    call cg_nsections_f(cg%unit, cg%base, zone, nsections, ier)
    if (ier .eq. ERROR) call cg_error_exit_f
!
    do s=1,nsections
      call cg_section_read_f(cg%unit, cg%base, zone, s, elsname, type, start, end, &
                             nbndry, parent_flag, ier)
      if (ier .eq. ERROR) call cg_error_exit_f
      esection(1,s) = start; esection(2,s) = end
      etype(s) = type
      call cg_npe_f(type,npe,ier)
      ne2vx(s) = npe
!
! read in element data now
      call cg_ElementDataSize_f(cg%unit, cg%base, zone, s, eldatasize, ier)
      if (ier .eq. ERROR) call cg_error_exit_f
!!write(*,*)'start/end=',s, trim(elsname), type, start,end, nbndry, parent_flag
      allocate(elements(1:eldatasize))
!!write(*,*) 'ne2vxmax,eldatasize,zone,ne2vx=',ne2vxmax,eldatasize,zone,ne2vx(s)
      call cg_elements_read_f(cg%unit, cg%base, zone, s, elements, null, ier)
      if (ier .eq. ERROR) call cg_error_exit_f
! move data into place
!!write(*,*) 'e2vx load begin...',nelem
      k=0
      do i=start,end
        do j=1,ne2vx(s)
          k=k+1
          e2vx(j,i)=elements(k)
        enddo
      enddo
!!write(*,*) 'e2vx stored for section=',s
      deallocate(elements)
    enddo
!
    if (ALLOCATED(elements)) deallocate(elements)

    end subroutine cgns_element_info
!
    subroutine cgns_cell_data(cg,x,y,z,nvx)
!
    implicit none
    integer :: nvx
    real :: x(nvx),y(nvx),z(nvx)
    type (cgnsdb) :: cg
!
    integer :: ncoords,coord
    integer :: nzones,zone,zonetype,datatype
    integer :: rmin,rmax
    character (len=32) :: coordname,zonename,name
!
! A cell is equal to a zone in the CGNS database
!
! Zone
    zone=1
!
!--- zone attribute:  GOTO zone node
    call cg_goto_f(cg%unit, cg%base, ier, 'Zone_t', zone, 'end')
    if (ier .eq. ERROR) call cg_error_exit_f
!
!--- zone coordinate attribute:  GOTO GridCoordinates_t node
    call cg_goto_f(cg%unit, cg%base, ier, 'Zone_t', zone,  &
                  'GridCoordinates_t', 1, 'end')
    if (ier .eq. ERROR) call cg_error_exit_f
!
!--- read coordinates using coordinate arrays' specific functions:
    call cg_ncoords_f(cg%unit, cg%base, zone, ncoords, ier)
    if (ier .eq. ERROR) call cg_error_exit_f
!
    rmin=1;rmax=nvx

    do coord=1, ncoords
      call cg_coord_info_f(cg%unit, cg%base, zone, coord, datatype, &
                                 coordname, ier)
      if (ier .eq. ERROR) call cg_error_exit_f
!
      if (trim(coordname)=='CoordinateX') then
        call cg_coord_read_f(cg%unit, cg%base, zone, coordname, &
                  RealDouble, rmin, rmax, x, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
      elseif (trim(coordname)=='CoordinateY') then
        call cg_coord_read_f(cg%unit, cg%base, zone, coordname, &
                    RealDouble, rmin, rmax, y, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
      elseif (trim(coordname)=='CoordinateZ') then
        call cg_coord_read_f(cg%unit, cg%base, zone, coordname, &
                     RealDouble, rmin, rmax, z, ier)
        if (ier .eq. ERROR) call cg_error_exit_f
      else
        write(*,*) 'Wrong coordinate input in cgns_cell_data'
        stop
      endif
    enddo

    end subroutine cgns_cell_data
!

end module
