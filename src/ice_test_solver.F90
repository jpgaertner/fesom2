MODULE ice_test_solver
  USE MOD_ICE
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_MESH
  USE o_PARAM
  USE g_comm_auto
  IMPLICIT NONE
  private
  public :: solve_test

contains

  subroutine solve_test(ice, partit, mesh)
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in), target :: mesh

    integer :: i, nod, edge
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, u_ice_nod, v_ice_nod
    real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
    real(kind=WP), dimension(:), pointer  :: u_w, v_w
    integer, pointer :: ice_vplace
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    ice_vplace => ice%ice_vplace  
    u_ice => ice%uice(:)
    v_ice => ice%vice(:)
    u_ice_nod => ice%uice_nod(:)
    v_ice_nod => ice%vice_nod(:)
    stress_atmice_x => ice%stress_atmice_x(:)
    stress_atmice_y => ice%stress_atmice_y(:)
    u_w => ice%srfoce_u(:)
    v_w => ice%srfoce_v(:)


    u_ice = 0.0_WP
    v_ice = 0.0_WP

    ! test if the exchange routines work
    if (ice_vplace == 0) then
      do i = 1, myDim_nod2D
        u_ice(i) = 1.0_WP
        v_ice(i) = 1.0_WP
      end do

      print *, 'u_ice nodes = ', int(u_ice(1)), int(u_ice(myDim_nod2D)), int(u_ice(myDim_nod2D+1)), int(u_ice(myDim_nod2D+eDim_nod2D))
      call exchange_nod(u_ice, partit)
      print *, 'u_ice nodes = ', int(u_ice(1)), int(u_ice(myDim_nod2D)), int(u_ice(myDim_nod2D+1)), int(u_ice(myDim_nod2D+eDim_nod2D))

    else if (ice_vplace == 1) then
      do i = 1, myDim_edge2D
        u_ice(i) = 1.0_WP
        v_ice(i) = 1.0_WP
      end do

      print *, 'u_ice edges = ', int(u_ice(1)), int(u_ice(myDim_edge2D)), int(u_ice(myDim_edge2D+1)), int(u_ice(myDim_edge2D+eDim_edge2D))
      call exchange_edge2D(u_ice, partit)
      print *, 'u_ice edges = ', int(u_ice(1)), int(u_ice(myDim_edge2D)), int(u_ice(myDim_edge2D+1)), int(u_ice(myDim_edge2D+eDim_edge2D))
    end if


    ! interpolate velocities from edges to nodes if necessary
    if (ice_vplace == 0) then
      u_ice_nod = u_ice
      v_ice_nod = v_ice
    else if (ice_vplace == 1) then
      u_ice_nod = 0
      v_ice_nod = 0

      do edge = 1, myDim_edge2D
        ! sum the adjacent edge velocities at the nodes
        u_ice_nod(edges(1,edge)) = u_ice_nod(edges(1,edge)) + u_ice(edge)
        u_ice_nod(edges(2,edge)) = u_ice_nod(edges(2,edge)) + u_ice(edge)
        v_ice_nod(edges(1,edge)) = v_ice_nod(edges(1,edge)) + v_ice(edge)
        v_ice_nod(edges(2,edge)) = v_ice_nod(edges(2,edge)) + v_ice(edge)
      end do

      do nod = 1, myDim_nod2D
        ! and divide by the number of adjacent edges
        u_ice_nod(nod) = u_ice_nod(nod) / mesh%nn_num(nod)
        v_ice_nod(nod) = v_ice_nod(nod) / mesh%nn_num(nod)
      end do

    end if

    !print *, 'myDim_nod2D = ', myDim_nod2D
    !print *, 'myDim_edge2D = ', myDim_edge2D
    !print *, 'myDim_elem2D = ', myDim_elem2D


  end subroutine solve_test


END MODULE ice_test_solver
