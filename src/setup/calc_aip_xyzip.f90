!
!  calc_xyzip_uns.f90
!
module mod_calc_aip_xyzip
contains
!
  subroutine calc_aip_xyzip_uns(x,y,z,xyzip,aip,e2vx,etype, &
                                 nelem,nbndry,nvx,ne2vxmax,f2e,nsec,esec,nf)
!
      use mod_util
      implicit none
!
      integer :: nelem,nvx,ne2vxmax,nsec,nf,nbndry
      real :: x(nvx),y(nvx),z(nvx),xyzip(3,nf),aip(3,nf)
      integer :: e2vx(ne2vxmax,nelem),etype(nsec),f2e(nf),esec(2,nsec)
!
      real :: r1(3),r2(3),r3(3),r4(3),dr1(3),dr2(3)
      real, dimension(3,1) :: subcntr,areavec,centeroid
      real, dimension(3,3) :: sumcntr
      real, dimension(3) :: tmp
      integer :: e,s,t,fl,fg,lst(ne2vxmax),nl,i
!
! xyzip calculation

    do fg=1,nf
        call get_idx(f2e(fg),0,e,fl)
        do s=1, nsec
            if(e>=esec(1,s) .and. e<=esec(2,s)) exit
        enddo
        t=etype(s)
        call face_vxlist(e2vx(:,e),t,fl,lst,nl)
        ! Process triangles
        if (nl==3) then
            r1(1) = x(lst(1)); r1(2) = y(lst(1)); r1(3)=z(lst(1))
            r2(1) = x(lst(2)); r2(2) = y(lst(2)); r2(3)=z(lst(2))
            r3(1) = x(lst(3)); r3(2) = y(lst(3)); r3(3)=z(lst(3))
!
            subcntr(:,1) = (r1+r2+r3)/(3.0)
            dr1=r2-r1;dr2=r3-r1
            call vec_a_bcrossc(areavec,dr1,dr2,1)
            areavec=(0.5)*areavec
            aip(:,fg)=areavec(:,1)
            sumcntr=matmul(subcntr,transpose(areavec))
        elseif (nl==4) then
            r1(1) = x(lst(1)); r1(2) = y(lst(1)); r1(3)=z(lst(1))
            r2(1) = x(lst(2)); r2(2) = y(lst(2)); r2(3)=z(lst(2))
            r3(1) = x(lst(3)); r3(2) = y(lst(3)); r3(3)=z(lst(3))
            r4(1) = x(lst(4)); r4(2) = y(lst(4)); r4(3)=z(lst(4))
!
            dr1=r2-r1;dr2=r3-r1
            call vec_a_bcrossc(areavec,dr1,dr2,1)
            areavec=(0.5)*areavec
            aip(:,fg)=areavec(:,1)
            subcntr(:,1) = (r1+r2+r3)/(3.0)
            sumcntr=matmul(subcntr,transpose(areavec))

            dr1=r3-r1;dr2=r4-r1
            call vec_a_bcrossc(areavec,dr1,dr2,1)
            areavec=(0.5)*areavec
            aip(:,fg)=aip(:,fg)+areavec(:,1)
            subcntr(:,1) = (r1+r3+r4)/(3.0)
            sumcntr=sumcntr+matmul(subcntr,transpose(areavec))

        else
            write(*,*) 'Unsupported face in calc_xyzip_uns ...'
            stop
        endif
        areavec(:,1)=aip(:,fg)
        centeroid=matmul(sumcntr,areavec)/dot_product(areavec(:,1),areavec(:,1))
        xyzip(:,fg) = centeroid(:,1)

    enddo
!
  end subroutine

  end module


