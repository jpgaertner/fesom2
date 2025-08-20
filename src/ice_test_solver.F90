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

  subroutine solve_test_old(ice, partit, mesh)
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in), target :: mesh
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: u_buff, v_buff
    integer       :: i, i1, i2
    real(kind=WP) :: expected_u, expected_v, diff_u, diff_v
    real(kind=WP) :: x1, x2, y1, y2, x_edge, y_edge
    !___________________________________________________________________________
    real(kind=WP), allocatable :: send_map(:)
    integer :: idx, edge, np, ierr, sender
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_buff => ice%nc%u_buff
    v_buff => ice%nc%v_buff

    u_buff = 0.0_WP
    v_buff = 0.0_WP

    do i = 1, myDim_edge2D
        i1 = edges(1, i)
        i2 = edges(2, i)
    
        x1 = coord_nod2D(1, i1)
        y1 = coord_nod2D(2, i1)
        x2 = coord_nod2D(1, i2)
        y2 = coord_nod2D(2, i2)
    
        x_edge = 0.5_wp * (x1 + x2)
        y_edge = 0.5_wp * (y1 + y2)
    
        !u_ice(i) = sin(x_edge) + cos(y_edge)
        !v_ice(i) = sin(x_edge) - cos(y_edge)

        ! use this one to catch errors in the orientation of edges
        u_buff(i) = x2 - x1
        v_buff(i) = y2 - y1
    end do

    call exchange_edge2D(u_buff, partit)
    call exchange_edge2D(v_buff, partit)

    do i = myDim_edge2D+1, myDim_edge2D+eDim_edge2D
        i1 = edges(1, i)
        i2 = edges(2, i)

        x1 = coord_nod2D(1, i1)
        y1 = coord_nod2D(2, i1)
        x2 = coord_nod2D(1, i2)
        y2 = coord_nod2D(2, i2)

        x_edge = 0.5_wp * (x1 + x2)
        y_edge = 0.5_wp * (y1 + y2)
    
        !expected_u = sin(x_edge) + cos(y_edge)
        !expected_v = sin(x_edge) - cos(y_edge)

        expected_u = x2 - x1
        expected_v = y2 - y1

        diff_u = abs(u_buff(i) - expected_u)
        diff_v = abs(v_buff(i) - expected_v)

        if (diff_u > 1e-10_WP .or. diff_v > 1e-10_WP) then
            print *, 'Halo mismatch at edge', i, diff_u, diff_v
        end if
    end do

  end subroutine solve_test_old

  subroutine solve_test(ice, partit, mesh)
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in), target :: mesh
    !___________________________________________________________________________
    real(kind=WP), pointer :: u_buff(:), v_buff(:)
    integer :: i, j, pe, edge_ghost, ed_global, edge_ghost_global
    real(kind=WP) :: expected_u, expected_v, diff_u, diff_v
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_buff => ice%nc%u_buff
    v_buff => ice%nc%v_buff
    
    ! initialize local edges with a unique pattern per PE/edge
    do i = 1, myDim_edge2D
     ed_global = myList_edge2D(i)
     u_buff(i) = 1000*mype + ed_global
     v_buff(i) = 2000*mype + ed_global
    end do

    ! fill ghost edges
    call exchange_edge2D(u_buff, partit)
    call exchange_edge2D(v_buff, partit)
    
    ! loop over received ghost edges using com_edge2D
    do pe = 1, com_edge2D%rPEnum
     do j = com_edge2D%rptr(pe), com_edge2D%rptr(pe+1)-1
        edge_ghost = com_edge2D%rlist(j)
        edge_ghost_global = myList_edge2D(edge_ghost)
    
        ! expected values come from the actual sending PE
        expected_u = 1000*com_edge2D%rPE(pe) + edge_ghost_global
        expected_v = 2000*com_edge2D%rPE(pe) + edge_ghost_global
    
        diff_u = abs(u_buff(edge_ghost) - expected_u)
        diff_v = abs(v_buff(edge_ghost) - expected_v)
    
        if (diff_u > 1e-10_WP .or. diff_v > 1e-10_WP) then
           print *, 'Halo mismatch at edge', edge_ghost_global, 'from PE', com_edge2D%rPE(pe), &
                    'diff_u=', diff_u, 'diff_v=', diff_v
        end if
     end do
    end do

end subroutine solve_test


END MODULE ice_test_solver
