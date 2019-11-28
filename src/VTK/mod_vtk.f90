module mod_vtk

#ifdef Catalyst
    INTERFACE
      SUBROUTINE cfd_catalyst(command, phic, grad, u, v, w, temp, px, py, pz, e2vx, npnts, ncvs, esec, etype, nsec, ne2vx_max, &
								 iout, max_it, dt) BIND(C)
        USE, INTRINSIC :: ISO_C_BINDING
        !character :: command
        integer :: command
        REAL  :: phic, grad, u, v, w, temp
        real  :: px,py,pz,dt
        integer  :: npnts,ncvs, es,ee,l, e2vx, esec, etype, nsec, ne2vx_max, iout, max_it

      END SUBROUTINE cfd_catalyst
    END INTERFACE
#endif

#ifdef VTK
    INTERFACE
      SUBROUTINE c_vtk_writer( phic, grad, u, v, w, temp, px, py, pz, e2vx, npnts, ncvs, esec, etype, nsec, ne2vx_max, it) BIND(C)
        USE, INTRINSIC :: ISO_C_BINDING
!        character :: filename
        REAL  :: phic, grad, u, v, w, temp
        real  :: px,py,pz
        integer  :: npnts,ncvs, es,ee,l, e2vx, esec, etype, nsec, ne2vx_max, it

      END SUBROUTINE c_vtk_writer
    END INTERFACE
#endif
contains
  subroutine vtk_data(eqn,eqnt,geom,iout)
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
       type(uvwp_t), target :: eqn
       type(energy_t), target :: eqnt
       integer :: iout, file_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VTK
        call c_vtk_writer( eqn%p(1), eqn%gpc(1), eqn%u(1), eqn%v(1), eqn%w(1), eqnt%t(1), geom%x(1), geom%y(1), geom%z(1), geom%mg%e2vx(1),&
						   geom%nvx,geom%ne, geom%mg%esec(1,1), geom%mg%etype(1),geom%mg%nsec,geom%mg%ne2vx_max, iout)
#endif
!        print*, "vtk write finished!"
  end subroutine vtk_data

end module




