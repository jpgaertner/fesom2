!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    


module testcase
  contains
  subroutine exchange_test(partit)
    use mod_partit
    use g_comm_auto
    implicit none
    type(t_partit), intent(inout), target :: partit
    real(kind=WP) :: uIce(partit%myDim_edge2D+partit%eDim_edge2D)
    real(kind=WP) :: vIce(partit%myDim_nod2D+partit%eDim_nod2D)
    integer :: i
#include "associate_part_def.h"
#include "associate_part_ass.h"

    uIce = 0
    do i=1,myDim_edge2D
      uIce(i) = 1
    end do

    vIce = 0
    do i=1,myDim_nod2D
      vIce(i) = 1
    end do

    print *, 'sometest', com_nod2D%rPEnum, com_edge2D%rPEnum

    print *, 'edge_init', int(uIce(1)), int(uIce(myDim_edge2D)), int(uIce(myDim_edge2D+1)), int(uIce(myDim_edge2D+eDim_edge2D))
    call exchange_edge2D(uIce, partit)
    print *, 'edge_exch', int(uIce(1)), int(uIce(myDim_edge2D)), int(uIce(myDim_edge2D+1)), int(uIce(myDim_edge2D+eDim_edge2D))

    print *, 'nod_init', int(vIce(1)), int(vIce(myDim_nod2D)), int(vIce(myDim_nod2D+1)), int(vIce(myDim_nod2D+eDim_nod2D))
    call exchange_nod(vIce, partit)
    print *, 'nod_exch', int(vIce(1)), int(vIce(myDim_nod2D)), int(vIce(myDim_nod2D+1)), int(vIce(myDim_nod2D+eDim_nod2D))



  end subroutine

end module

program main
  use fesom_module
  use fesom_main_storage_module
  use testcase

  integer nsteps

  call fesom_init(nsteps)
  call exchange_test(f%partit)
  call fesom_runloop(nsteps)
  call fesom_finalize

end program