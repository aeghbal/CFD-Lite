module mod_vtu_output
  contains
!
! write_vtubin
!
       subroutine write_vtubin(eqn,geom,projPath,iout)
       use mod_eqn_setup
       use mod_cell
       use mod_util
       implicit none

       character(len=180) :: vtu_filename,projPath,outputPath
       type(geometry_t) :: geom
       class(equation_t) :: eqn

       integer :: iout,file_unit
       integer :: alength,anlength,i,j,k,l,m,n,ifile,ivar,type_no
       character (len=8) :: ctype,string,cit
       character (len=12) :: cformat,cnewline
       character (len=180) :: celli,inti,cval,string1,string2,string3,filename,filename0,casename,stroffset,oper,cappen,cexpbc

       character (len=180) :: vari,cbcfamname,cfolder,var
       integer (kind=8), pointer, dimension(:) :: nbytes2,nbytes
       character (len=1) :: lfed

! make output directory
      write (string, '(I8)') iout
      outputPath=trim(projPath)//'VTK/output-'//trim(adjustl(string))
      call check_casedir(outputPath)
      vtu_filename=trim(projPath)//'VTK/output-'//trim(adjustl(string))//'/'//trim(eqn%name)//'.vtu'

      file_unit=111
      open ( unit=file_unit, status='unknown',access='stream', form='unformatted' , file=trim(vtu_filename))
      !open ( unit=file_unit, status='unknown', form='binary', file=trim(vtu_filename))
! The VTU file for cell data

      oper = 'pre-cell'
      call vtu_data(oper,eqn,geom,file_unit,nbytes,nbytes2)

      oper = 'write-cldata'
      call vtu_data(oper,eqn,geom,file_unit,nbytes,nbytes2)

      close(file_unit)


      deallocate(nbytes,nbytes2)

    end subroutine write_vtubin
!
! vtu_data
!
       subroutine vtu_data(oper,eqn,geom,file_unit,nbytes,nbytes2)
       use mod_cell
       use mod_util
       use mod_eqn_setup
       use mod_uvwp
       use mod_energy
       use mod_scalar
       use mod_vfr
       use mod_mfr
       implicit none

       type(geometry_t) :: geom
       class(equation_t), target :: eqn
       integer :: file_unit


       integer :: icvs,jcvs,kcvs,type_no,fg,v,idx,fg_sgn,e
       integer :: alength,anlength,i,j,k,l,m,n,ifile,ivar,icl,ibase,ii,jj,kk
       character (len=8) :: ctype,string
       character (len=12) :: cformat,cnewline
       character (len=32) :: CLfamilyname
       character (len=180) :: celli,inti,cval,string1,string2,string3,filename,filename0,casename,stroffset,suffix
       character (len=*) :: oper

       integer :: itemp,iqp,ipass,sz,ncvs,noffset,npnts
       integer :: nvarlist,ntrvarbclist,nbcfam,tint,length,cnt
       character (len=180) :: vari,cbcfamname,checkvar
       character (len=32), allocatable :: varname(:)
       character (len=180), allocatable :: varlabel(:),vartype(:), &
                                          varmemory(:),vtkarraytype(:)
       character (len=60), allocatable :: varmodel(:)
       character (len=1) :: lfed
       real , allocatable :: vardimension(:)
       real :: tfloat
       !integer (kind=8) :: nbytes(*),nbytes2(*)
       integer (kind=8), pointer, dimension(:) :: nbytes2,nbytes
       integer (kind=8) :: tint8

       integer, pointer :: ipts(:),gf2vx(:),gf2vx_idx(:)
       type :: AoD_t! array of double
          real, pointer :: p(:)
       end type
       type(AoD_t), allocatable :: phi(:)


    noffset = 0
    lfed = char(10)
    cformat = 'appended'

! Get output variable list
    select type(eqn)

      type is(uvwp_t)
        nvarlist=7
        allocate(phi(nvarlist))
        allocate (varname(nvarlist))
        allocate (vardimension(nvarlist))
        allocate(vtkarraytype(nvarlist))
        phi(1)%p=>eqn%u;varname(1)='u';vardimension(1)=1
        phi(2)%p=>eqn%v;varname(2)='v';vardimension(2)=1
        phi(3)%p=>eqn%w;varname(3)='w';vardimension(3)=1
        phi(4)%p=>eqn%p;varname(4)='p';vardimension(4)=1
        phi(5)%p=>eqn%gpc;varname(5)='gpc';vardimension(5)=3
        do e=1,geom%ne
          eqn%dc(e)=0.
          do idx=geom%ef2nb_idx(e),geom%ef2nb_idx(e+1)-1
            fg=geom%ef2nb(idx,2)
            fg_sgn=fg/abs(fg)
            fg=abs(fg)
            eqn%dc(e)=eqn%dc(e)+eqn%mip(fg)*fg_sgn
          enddo
        enddo
        phi(6)%p=>eqn%dc;varname(6)='mip';vardimension(6)=1
!       !!!!!!!!SUBD!!!!!
!        do c=1,nsubd
!          do e=1,subdomain(c)%ne
!            idx=geom%mg%g2gf(1)%idx(c)+e-1
!            g=geom%mg%g2gf(1)%p(idx)
!            phi(g) = subdomain(c)%phic(e)
!          enddo
!        enddo
      type is(energy_t)
        nvarlist=4
        allocate(phi(nvarlist))
        allocate (varname(nvarlist))
        allocate (vardimension(nvarlist))
        allocate(vtkarraytype(nvarlist))
        phi(1)%p=>eqn%phi;varname(1)='enthalpy';vardimension(1)=1
        phi(2)%p=>eqn%grad;varname(2)='grad_enthalpy';vardimension(2)=3
        phi(3)%p=>eqn%t;varname(3)='temperature';vardimension(3)=1
        phi(4)%p=>eqn%gt;varname(4)='grad_temperature';vardimension(4)=3
      class default
        nvarlist=2
        allocate(phi(nvarlist))
        allocate (varname(nvarlist))
        allocate (vardimension(nvarlist))
        allocate(vtkarraytype(nvarlist))
        phi(1)%p=>eqn%phi;varname(1)='phi';vardimension(1)=1
        phi(2)%p=>eqn%grad;varname(2)='grad';vardimension(2)=3
    end select
    vtkarraytype(1:nvarlist)='Float64'

    if(trim(oper)=='pre-cell') then
       allocate (nbytes(nvarlist))
       allocate (nbytes2(6))
       nbytes = 0 ; nbytes2 = 0
    endif

    npnts = geom%nvx
    ncvs = geom%ne

    if(trim(oper)=='pre-cell') then
      write(file_unit) '<?xml version="1.0"?>'//lfed
      write(file_unit) '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lfed
      write(file_unit) ' <UnstructuredGrid>'//lfed
      write (string1, '(I12)') npnts
      write (string2, '(I12)') ncvs
      write(file_unit) '  <Piece NumberOfPoints="'//trim(adjustl(string1))//'" NumberOfCells="'//trim(adjustl(string2))//'">'//lfed
      write(file_unit) '   <CellData Scalars="scalars">'//lfed

      do ivar = 1, nvarlist

        checkvar = varname(ivar)

        write (string3, '(I1)') int(abs(vardimension(ivar)))
        write (stroffset, '(I12)') noffset

        write(file_unit) '    <DataArray type="'//trim(vtkarraytype(ivar))//'" Name="'//trim(varname(ivar))// &
                        '" NumberOfComponents="'//trim(string3)//'" format="'//trim(cformat)//'" offset="'//trim(stroffset)//'" />'//lfed
!
        cnt = ncvs

! Variables in single memory location
        nbytes(ivar) = nbytes(ivar) + cnt*int(abs(vardimension(ivar)))*sizeof(tfloat)

        noffset = noffset + sizeof(tint8) + nbytes(ivar)

      enddo

      write(file_unit) '   </CellData>'//lfed
      write (stroffset, '(I12)') noffset

! Mesh point coordinate Data
      write(file_unit) '   <Points>'//lfed
      write(file_unit) '    <DataArray type="Float64" NumberOfComponents="3" format="'//trim(cformat)//'" offset="'//trim(stroffset)//'" />'//lfed

      anlength = geom%nvx
      nbytes2(1) = nbytes2(1) + 3*anlength*sizeof(tfloat)
! Point connectivity of Cells
      noffset = noffset + sizeof(tint8) + nbytes2(1)
      write (stroffset, '(I12)') noffset

      write(file_unit) '   </Points>'//lfed
      write(file_unit) '   <Cells>'//lfed
      write(file_unit) '    <DataArray type="Int32" Name="connectivity" format="'//trim(cformat)//'" offset="'//trim(stroffset)//'" />'//lfed

      cnt=0
      do n=1,geom%mg%nsec
        if(ElementTypeDim(geom%mg%etype(n))==2) cycle!boundary
        cnt=cnt+(geom%mg%esec(2,n)-geom%mg%esec(1,n)+1)*element_nvx(geom%mg%etype(n))
      enddo
      nbytes2(2) = nbytes2(2) + cnt*sizeof(tint)

      noffset = noffset + sizeof(tint8) + nbytes2(2)
      write (stroffset, '(I12)') noffset

      write(file_unit) '    <DataArray type="Int32" Name="offsets" format="'//trim(cformat)//'" offset="'//trim(stroffset)//'" />'//lfed

      cnt=0
      do n=1,geom%mg%nsec
        if(ElementTypeDim(geom%mg%etype(n))==2) cycle!boundary
        cnt=cnt+geom%mg%esec(2,n)-geom%mg%esec(1,n)+1
      enddo

      nbytes2(3) = nbytes2(3) + cnt*1*sizeof(tint)

      noffset = noffset + sizeof(tint8) + nbytes2(3)
      write (stroffset, '(I12)') noffset

      write(file_unit) '    <DataArray type="Int32" Name="types" format="'//trim(cformat)//'" offset="'//trim(stroffset)//'" />'//lfed
!
      cnt=geom%ne
      nbytes2(4) = nbytes2(4) + cnt*1*sizeof(tint)
!
      noffset = noffset + sizeof(tint8) + nbytes2(4)

      write(file_unit) '   </Cells>'//lfed
      write(file_unit) '  </Piece>'//lfed

      write(file_unit) ' </UnstructuredGrid>'//lfed
      write(file_unit) ' <AppendedData encoding="raw">'//lfed
      write(file_unit) '_'

    else

       do ivar = 1, nvarlist

         checkvar = varname(ivar)

         write(file_unit) nbytes(ivar)
!
         cnt = ncvs

! Variables in single memory location
         cval=varname(ivar)
         k = vardimension(ivar)

         do n=1,geom%mg%nsec
            if(ElementTypeDim(geom%mg%etype(n))==2) cycle!boundary
            do m=geom%mg%esec(1,n),geom%mg%esec(2,n)
              write(file_unit) phi(ivar)%p((m-1)*k+1:m*k)
            enddo
         enddo
       enddo

       write(file_unit) nbytes2(1)

       do j = 1, geom%nvx
         write(file_unit) geom%x(j),geom%y(j),geom%z(j)
       enddo

       write(file_unit) nbytes2(2)

      cnt=0
      do n=1,geom%mg%nsec
        if(ElementTypeDim(geom%mg%etype(n))==2) cycle!boundary
        l=element_nvx(geom%mg%etype(n))
        do m=geom%mg%esec(1,n),geom%mg%esec(2,n)
          v=(m-1)*geom%mg%ne2vx_max
          write(file_unit) geom%mg%e2vx(v+1:v+l)-1! -1 is for the offset. the vertex indeces start from zero
        enddo
      enddo

      l = 0
      write(file_unit) nbytes2(3)
      do n=1,geom%mg%nsec
        if(ElementTypeDim(geom%mg%etype(n))==2) cycle!boundary
        do m=geom%mg%esec(1,n),geom%mg%esec(2,n)
          l=l+element_nvx(geom%mg%etype(n))
          write(file_unit) l
        enddo
      enddo

      write(file_unit) nbytes2(4)

      do j = 1, geom%mg%nsec
        if(ElementTypeDim(geom%mg%etype(j))==2) cycle!boundary
        select case(geom%mg%etype(j))
          case(10)
            type_no=10!tetrahedron
          case(17)
            type_no=12!hexahedron
          case(12)
            type_no=14!pyramid
          case(14)
            type_no=13 !prism or wedge
          case(23)
            type_no=42 ! polyhedron
          case default
            print *,'Unknown uns element type...' !other formats
            stop
        end select
        do k = geom%mg%esec(1,j),geom%mg%esec(2,j)
          write(file_unit) type_no
        enddo
      enddo

      write(file_unit) lfed//' </AppendedData>'//lfed
      write(file_unit) '</VTKFile>'//lfed
    endif
!
    if (ALLOCATED(varlabel)) deallocate(varlabel)
    if (ALLOCATED(vartype)) deallocate(vartype)
    if (ALLOCATED(varmemory)) deallocate(varmemory)
    if (ALLOCATED(varname)) deallocate(varname)
    if (ALLOCATED(vtkarraytype)) deallocate(vtkarraytype)
    if (ALLOCATED(varmodel)) deallocate(varmodel)
    if (ALLOCATED(vardimension)) deallocate(vardimension)
    do ivar=1,nvarlist
      nullify(phi(ivar)%p)
    enddo
    if (allocated(phi)) deallocate(phi)

  end subroutine vtu_data

end module
