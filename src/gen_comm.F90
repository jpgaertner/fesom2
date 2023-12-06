! Cell-vertex finite-volume version
! Contains: Routines that support parallelization
! set_par_support_ini run in the initialization phase.
! The communication rules are saved. 
! set_par_support in the main phase just allocates memory for buffer 
! arrays, the rest is read together with mesh from saved files.
!=======================================================================
subroutine communication_nodn(partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                  :: n, np, prank, el, r_count, s_count, q, i, j, nod, k, l
  integer                  :: num_send(0:partit%npes-1), num_recv(0:partit%npes-1), nd_count
  integer, allocatable     :: recv_from_pe(:), send_to_pes(:,:)
  logical                  :: max_laendereck_too_small=.false.
  integer                  :: IERR
#include "associate_part_def.h"
#include "associate_mesh_ini.h"
#include "associate_part_ass.h" !part only
  ! Assume we have 2D partitioning vector in part. Find communication rules
  ! Reduce allocation: find all neighboring PE
  nd_count = count(part(1:nod2d) == mype)
! write(*,*) nod2d
! write(*,*) MAX_LAENDERECK
! write(*,*) nd_count
! write(*,*) allocated(partit%myList_nod2D)
! write(*,*) partit%mype
  allocate(recv_from_pe(nod2d), send_to_pes(MAX_LAENDERECK,nd_count), &
           partit%myList_nod2D(nd_count), STAT=IERR)
           !# ??? what is ierr? it is not assigned to something in the
           ! associate files 
           !# ??? myList_nod2D does not contain halo nodes?
           ! how can it be used like myList_nod2D(1:myDim_nod2D+eDim_nod2D)?
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_nodn'
     stop
  endif

  myList_nod2D=>partit%myList_nod2D
  nd_count = 0
  do n=1,nod2D
     ! Checks if element el has nodes that belong to different partitions
     if (part(n) == mype) then
        nd_count = nd_count+1
        myList_nod2D(nd_count)=n
     endif
  end do
  myDim_nod2D=nd_count

  num_send(0:npes-1) = 0
  num_recv(0:npes-1) = 0
  recv_from_pe(1:nod2d) = -1
  send_to_pes(1:MAX_LAENDERECK,1:nd_count) = -1

  ! For the local nodes, run through the patch and collect all nodes in the patch
  ! (It would be simpler to re-use the adjacency matrix, but it is not a global
  !  variable... and I am lazy and want to keep the data structure simple)

  do l=1,nd_count
     n = myList_nod2D(l)
     do i = 1, nod_in_elem2D_num(n)
      !# ??? nod_in_elem2D_num(n) is the number of elements that contain the node n
        ! Over all elements el that the node n is part of
        el = nod_in_elem2D(i,n)
        !# ??? nod_in_elem2D(:,n) is the elements that the node n is part of
        ! el is now the i-th element that uses node n

        ! Checks, if elements are quads or triangles
        q = 4    ! quads as default
        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

        do j = 1, q
           ! Over all nodes in every element el
           nod = elem2D_nodes(j,el)
           
           ! Checks, if node j is not in another partitionen
           if (part(nod) /= mype) then
              ! Checks, if not already considered to be received from this
              ! node from partition part(nod)
              if (recv_from_pe(nod) == -1) then  ! nod already collected to be received?
                 ! We have to receive this node from part(nod)
                 ! Add plus one node to the total number of
                 num_recv(part(nod)) = num_recv(part(nod)) + 1
                 ! ???
                 recv_from_pe(nod) = part(nod)   ! recv_from_pe(recv_count) = nod  ! no new information, just handy
              endif
              ! Checks, if all possible connected partition
              ! And we have to send n to part(nod). Do we know this already?
              do k=1,MAX_LAENDERECK    !???
               !# ???
                 if (send_to_pes(k,l) == part(nod)) then
                    exit  ! already collected
                 elseif (send_to_pes(k,l) == -1) then
                    send_to_pes(k,l) = part(nod)
                    num_send(part(nod)) = num_send(part(nod)) + 1
                    exit
                 elseif (k== MAX_LAENDERECK) then
                    max_laendereck_too_small = .true.  ! Problem
                 endif
              enddo
           endif
        enddo
     enddo
  enddo

  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Now, build the send and recv communication data structure
  ! To how many PE needs to be send and which PEs
!!$  com_nod2D%rPEnum = count(num_recv(0:npes-1) > 0)
!!$  com_nod2D%sPEnum = count(num_send(0:npes-1) > 0)

!!$  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
!!$       com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
!!$     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
!!$     stop
!!$  endif
!!$  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
!!$  allocate(com_nod2D%sPE(com_nod2D%sPEnum))

  r_count = 0
  s_count = 0
  com_nod2D%rptr(1) = 1
  com_nod2D%sptr(1) = 1

  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_nod2D%rPE(r_count) = np
        com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_nod2D%sPE(s_count) = np
        com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
     end if
  enddo
  com_nod2D%rPEnum = r_count
  com_nod2D%sPEnum = s_count
  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
       com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Counts the number of node for each partition PE mype has to send/receive
  ! In ascending order of PE number
!!$  r_count = 0
!!$  s_count = 0
!!$  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
!!$  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))

!!$  com_nod2D%rptr(1) = 1
!!$  com_nod2D%sptr(1) = 1
!!$
!!$  do r_count = 1, com_nod2D%rPEnum
!!$     np = com_nod2D%rPE(r_count)
!!$     com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
!!$  enddo
!!$  do s_count = 1, com_nod2D%sPEnum
!!$     np = com_nod2D%sPE(s_count)
!!$     com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
!!$  enddo

  ! Lists themselves

  r_count = 0
  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
  allocate(partit%com_nod2D%rlist(eDim_nod2D), &
           partit%com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1), STAT=IERR) 
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_nodn'
     stop
  endif

  do np = 1,com_nod2D%rPEnum
     prank = com_nod2D%rPE(np)
     do n = 1, nod2D
        if (recv_from_pe(n) == prank) then
           r_count = r_count+1
           com_nod2D%rlist(r_count) = n
        end if
     end do
  end do

  s_count = 0
!!$  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np = 1,com_nod2D%sPEnum
     prank = com_nod2D%sPE(np)
     do l = 1, nd_count
        n = myList_nod2D(l)
        if(any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           com_nod2D%slist(s_count) = n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part.

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)

  deallocate(recv_from_pe, send_to_pes)
end subroutine communication_nodn

!==========================================================================
subroutine communication_elemn(partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer, allocatable     :: recv_from_pe(:), send_to_pes(:,:)
  logical                  :: max_laendereck_too_small=.false.
  integer                  :: n, k, ep, np, prank, el, nod
  integer                  :: p, q, j, elem, i, l, r_count, s_count, el_count
  integer                  :: num_send(0:partit%npes-1), num_recv(0:partit%npes-1)
  integer                  :: IERR
#include "associate_part_def.h"
#include "associate_mesh_ini.h"
#include "associate_part_ass.h" !part only
  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An elem is external to element n if neither of its nodes 
  ! belongs to PE, but it is among the neighbors. Element n belongs to PE if 
  ! any of its nodes does. 
  
  ! This routine takes into account 
  ! com_elem2D_full: all  neighbors  
  ! com_elem2D:      only those sharing an edge 


  !===========================================
  !  com_elem2D
  !===========================================
  com_elem2D     =>partit%com_elem2D
  com_elem2D_full=>partit%com_elem2D_full

  allocate(recv_from_pe(elem2D), STAT=IERR)
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_elemn'
     stop
  endif

  el_count = 0
  do el=1,elem2D
     ! Checks if element el has nodes that belong to different partitions
     if (any(part(elem2D_nodes(1:4,el)) == mype)) then
        el_count = el_count+1
        recv_from_pe(el_count) = el
     endif
  end do
  myDim_elem2D=el_count

  allocate(partit%myList_elem2D(el_count), send_to_pes(MAX_LAENDERECK,el_count), STAT=IERR)
  !# ??? this way myList_elem2D does not contain halo elements
  if (IERR /= 0) then
     write (*,*) 'Could not allocate arrays in communication_elemn'
     stop
  endif
  myList_elem2D=>partit%myList_elem2D

  myList_elem2D(1:el_count) = recv_from_pe(1:el_count)
  num_send(0:npes-1) = 0
  num_recv(0:npes-1) = 0
  recv_from_pe(1:elem2D) = -1
  send_to_pes(1:MAX_LAENDERECK,1:el_count) = -1
  
  ! For the local elements, collect all adjacent (sharing an edge) elements
  ! that belong to other PEs
  
  do l=1,el_count
     el = myList_elem2D(l)

     ! Checks if triangles or quads
     q = merge(3, 4, elem2D_nodes(1,el) == elem2D_nodes(4,el))

     do n = 1, q                       ! cycle through all nodes

        ! Neighboring element elem that shares an edge with element el
        elem = elem_neighbors(n,el)

        if (elem < 1) cycle  ! boundary, "ghost element"

        ! Check for elements to be received
        if (all(part(elem2D_nodes(1:4,elem)) /= mype) .and. recv_from_pe(elem)==-1) then
           ! elem to be received already collected?
           ! We have to receive elem from PE ep:
           ep = part(elem2D_nodes(1,elem))  ! PE of first node is "main" owner
           num_recv(ep) = num_recv(ep) + 1
           recv_from_pe(elem) = ep  
        endif

        ! Check for elements to be sent to
        ! And maybe, we have to send el to the owners of elem
        ! 1. Is partition mype the main owner of element el?
        if (part(elem2D_nodes(1,el)) == mype .and. any(part(elem2D_nodes(1:4,elem)) /= mype)) then

           ! 2. Who owns element elem and needs to get element el? We must check all nodes!
           p=merge(3, 4, elem2D_nodes(1,elem) == elem2D_nodes(4,elem))

           do i=1,p
              ep = part(elem2D_nodes(i,elem))

              ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
              if (any(part(elem2D_nodes(1:q,el)) == ep)) cycle

              ! 4. Ok, for the owner ep, check if sending el is already collected
              do k=1,MAX_LAENDERECK
                 if (send_to_pes(k,l) == ep) then
                    exit  ! already collected
                 elseif (send_to_pes(k,l) == -1) then
                    send_to_pes(k,l) = ep
                    num_send(ep) = num_send(ep) + 1
                    exit
                 elseif (k== MAX_LAENDERECK) then
                    max_laendereck_too_small = .true.  ! Problem
                 endif
              enddo
           enddo
        end if
     end do
  enddo

  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif
    
! Now, build the send and recv communication data structure
  r_count = 0
  s_count = 0
  com_elem2D%rptr(1) = 1
  com_elem2D%sptr(1) = 1

  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_elem2D%rPE(r_count) = np
        com_elem2D%rptr(r_count+1) =  com_elem2D%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_elem2D%sPE(s_count) = np
        com_elem2D%sptr(s_count+1) =  com_elem2D%sptr(s_count)+ num_send(np)
     end if
  enddo

  com_elem2D%rPEnum = r_count
  com_elem2D%sPEnum = s_count
  if (com_elem2D%rPEnum > MAX_NEIGHBOR_PARTITIONS .or.  &
       com_elem2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Lists themselves

  r_count = 0
  eDim_elem2D=com_elem2D%rptr(com_elem2D%rPEnum+1)-1   
  allocate(partit%com_elem2D%rlist(eDim_elem2D))
  do np = 1,com_elem2D%rPEnum
     prank = com_elem2D%rPE(np)
     do el = 1, elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           com_elem2D%rlist(r_count) = el
        end if
     end do
  end do
  
  s_count = 0
  allocate(partit%com_elem2D%slist(com_elem2D%sptr(com_elem2D%sPEnum+1)-1)) 
  do np = 1,com_elem2D%sPEnum
     prank = com_elem2D%sPE(np)
     do l = 1, el_count
        el = myList_elem2D(l)
        if( any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           com_elem2D%slist(s_count) = el
        end if
     end do
  end do
  
  !===========================================
  !  com_elem2D_full
  !===========================================

  ! The halo relations that are already determined can be kept.
  ! Just add the elements connected only via nodes.
  ! num_send(0:npes-1) = 0
  ! num_recv(0:npes-1) = 0
  ! recv_from_pe(1:elem2D) = -1
  ! send_to_pes(1:MAX_LAENDERECK,1:elem2D) = -1
  
  ! For the local elements, collect all adjacent (sharing a node) elements
  ! that belong to other PEs

  do l=1,el_count
     el = myList_elem2D(l)

        ! Checks if triangles or quads
        q = merge(3, 4, elem2D_nodes(1,el) == elem2D_nodes(4,el))

        do n = 1, q                       ! cycle through all nodes

           nod = elem2D_nodes(n,el)

           ! Loop over all elements that belong to node nod 
           do j = 1, nod_in_elem2D_num(nod)  ! and for each node, through its patch
              elem = nod_in_elem2D(j,nod)

              ! Check for elements to be received
              if (all(part(elem2D_nodes(1:4,elem)) /= mype) .and. recv_from_pe(elem)==-1) then
                 ! elem to be received already collected?
                 ! We have to receive elem from PE ep:
                 ep = part(elem2D_nodes(1,elem))  ! PE of first node is "main" owner
                 num_recv(ep) = num_recv(ep) + 1
                 recv_from_pe(elem) = ep  
              endif

              ! Check for elements to be sent to
              ! And maybe, we have to send el to the owners of elem
              ! This gets more complicated:
              ! 1. Is partition mype the main owner of element el?
              if (part(elem2D_nodes(1,el)) == mype .and. any(part(elem2D_nodes(1:4,elem)) /= mype)) then

                 ! 2. Who owns element elem and needs to get element el? We must check all nodes!
                 p=merge(3, 4, elem2D_nodes(1,elem) == elem2D_nodes(4,elem))

                 do i=1,p
                    ep = part(elem2D_nodes(i,elem))

                    ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
                    if (any(part(elem2D_nodes(1:q,el)) == ep)) cycle

                    ! 4. Ok, for the owner ep, check if sending el is already collected
                    do k=1,MAX_LAENDERECK
                       if (send_to_pes(k,l) == ep) then
                          exit  ! already collected
                       elseif (send_to_pes(k,l) == -1) then
                          send_to_pes(k,l) = ep
                          num_send(ep) = num_send(ep) + 1
                          exit
                       elseif (k== MAX_LAENDERECK) then
                          max_laendereck_too_small = .true.  ! Problem
                       endif
                    enddo
                 end do
              endif
           end do
        end do
  enddo
  
  if (max_laendereck_too_small) then
     print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
     stop
  endif

  ! Now, build the send and recv communication data structure

  r_count = 0
  s_count = 0
  com_elem2D_full%rptr(1) = 1
  com_elem2D_full%sptr(1) = 1

  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_elem2D_full%rPE(r_count) = np
        com_elem2D_full%rptr(r_count+1) =  com_elem2D_full%rptr(r_count)+ num_recv(np)
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_elem2D_full%sPE(s_count) = np
        com_elem2D_full%sptr(s_count+1) =  com_elem2D_full%sptr(s_count)+ num_send(np)
     end if
  enddo

  com_elem2D_full%rPEnum = r_count
  com_elem2D_full%sPEnum = s_count

  ! Lists themselves

  r_count = 0
  allocate(partit%com_elem2D_full%rlist(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1)) 
  do np = 1,com_elem2D_full%rPEnum
     prank = com_elem2D_full%rPE(np)
     do el = 1, elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           com_elem2D_full%rlist(r_count) = el
        end if
     end do
  end do

  s_count = 0
  allocate(com_elem2D_full%slist(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1)) 
  do np = 1,com_elem2D_full%sPEnum
     prank = com_elem2D_full%sPE(np)
     do l = 1, el_count
        el = myList_elem2D(l)
        if( any(send_to_pes(:,l) == prank)) then 
           s_count = s_count+1
           com_elem2D_full%slist(s_count) = el
        end if
     end do
  end do

  deallocate(recv_from_pe, send_to_pes)
end subroutine communication_elemn
!==========================================================================
subroutine communication_edgen(partit, mesh)
   use MOD_MESH
   use MOD_PARTIT
   use MOD_PARSUP
   implicit none
   type(t_mesh), intent(in), target        :: mesh
   type(t_partit), intent(inout), target   :: partit
   integer                 :: num_send(0:partit%npes-1), num_recv(0:partit%npes-1)
   integer, allocatable    :: recv_from_pe(:), send_to_pes(:,:)
   integer                 :: elem, k, k_ed, np, l, n
   integer                 :: elnodes(3), eledges(3)
   integer                 :: r_count, s_count, ed_count, prank
   logical                 :: max_laendereck_too_small=.false.
   ! ==========
   integer, pointer        :: rPEnum, sPEnum, rPE(:), sPE(:), rptr(:), sptr(:), rlist(:), slist(:)
#include "associate_part_def.h"
#include "associate_mesh_ini.h"
#include "associate_part_ass.h"
   ! initialize a pointer com_edge2D of type com_struct in associate_part_ass
   ! and assign it to partit%com_edge2D in associate_part_def
   ! com_edge2D of type com_struct is added to the type t_partit in mod_partit

   rPEnum  => com_edge2D%rPEnum
   sPEnum  => com_edge2D%sPEnum
   rPE     => com_edge2D%rPE
   sPE     => com_edge2D%sPE
   rptr    => com_edge2D%rptr
   sptr    => com_edge2D%sptr
   rlist   => com_edge2D%rlist
   slist   => com_edge2D%slist

   ! ==========
   ! the lists recv_from_pe/ send_to_pes contain the PEs an edge is received from/ send to
   ! ==========
   allocate(recv_from_pe(edge2D), send_to_pes(MAX_LAENDERECK,edge2D))
   num_recv(0:npes-1) = 0
   num_send(0:npes-1) = 0
   recv_from_pe(1:edge2D) = -1
   send_to_pes(1:MAX_LAENDERECK,1:edge2D) = -1

   do elem = 1, elem2D
       elnodes = elem2D_nodes(:,elem)
       eledges = elem_edges(:,elem)

       ! receive 
       do k = 1, 3
           ! check if node is in mype
           if (part(elnodes(k)) == mype) then
               ! only receive opposite edge if none of its nodes is in mype
               k_ed = modulo(k,3) + 1 ! this is the index of the edge opposite to the node k
               if (part(edges(1,eledges(k_ed))) /= mype .and. part(edges(2,eledges(k_ed))) /= mype) then
                   ! check if the edge is still not collected
                   if (recv_from_pe(eledges(k_ed)) == -1) then
                       np = part(edges(1,eledges(k_ed)))  ! PE of the first node of the edge
                       num_recv(np) = num_recv(np) + 1 ! number of received edges from PE np
                       recv_from_pe(eledges(k_ed)) = np   ! PE from which the edge is received
                   end if
               end if
           end if
       end do

       ! send
       do k = 1, 3
           np = part(elnodes(k))

           ! check if node is outside of mype
           if (np /= mype) then
               ! only send edge opposite of the node if it is not already in np
               k_ed = modulo(k,3) + 1 ! this is the index of the edge opposite to the node k
               if (part(edges(1,eledges(k_ed))) == np .or. part(edges(2,eledges(k_ed))) == np) cycle
               ! check if the edge is in mype
               ! only check for edges(1,..) so that the edge is not send from two PEs
               if (part(edges(1,eledges(k_ed))) == mype) then
                   do l = 1, MAX_LAENDERECK
                       ! check if the edge is already collected
                       if (send_to_pes(l,eledges(k_ed)) == np) then
                           exit
                       ! only add edge to send if it is not already collected
                       else if (send_to_pes(l,eledges(k_ed)) == -1) then
                           send_to_pes(l,eledges(k_ed)) = np  ! PE the edge is send to
                           num_send(np) = num_send(np) + 1 ! number of edges send to PE np
                           exit
                       else if (l == MAX_LAENDERECK) then
                           max_laendereck_too_small = .true.
                           ! if the loop is not exited before this, the number of neighbouring
                           ! PEs is too large
                       end if
                   end do
               end if
           end if
       end do

   end do

   if (max_laendereck_too_small) then
       print *,'Increase MAX_LAENDERECK in gen_modules_partitioning.F90 and recompile'
       stop
   endif

   ! ==========
   ! the number of PEs information is received from/ send to
   ! ==========
   rPEnum = count(num_recv(0:npes-1) > 0)
   sPEnum = count(num_send(0:npes-1) > 0)

   ! ==========
   ! the lists rPE/ sPE contain the PEs (numbered 0 to npes-1) information is received from/ send to
   ! ==========
   r_count = 0
   s_count = 0

   do np = 0, npes-1
       if(num_recv(np) /= 0) then  ! if something is received from the PE...
          r_count = r_count + 1
          rPE(r_count) = np        ! ... add it to the list rPE at the next index
       end if
       if(num_send(np) /= 0) then  ! if something is send to the PE...
          s_count = s_count + 1
          sPE(s_count) = np        ! ... add it to the list sPE at the next index
       end if
    enddo

    ! ==========
    ! the lists rptr/ sptr are used as pointers
    ! ==========
    r_count = 0
    s_count = 0

   rptr(1) = 1
   sptr(1) = 1

   do r_count = 1, rPEnum
       np = rPE(r_count)
       rptr(r_count+1) =  rptr(r_count) + num_recv(np)
   enddo
   do s_count = 1, sPEnum
       np = sPE(s_count)
       sptr(s_count+1) =  sptr(s_count) + num_send(np)
   enddo

   ! ==========
   ! the lists rlist/ slist contain the edges that are received/ send
   ! ==========
   allocate(com_edge2D%rlist(rptr(rPEnum+1)-1)) 
   allocate(com_edge2D%slist(sptr(sPEnum+1)-1))
   r_count = 0
   s_count = 0

   do np = 1,rPEnum
      prank = rPE(np)

      do n = 1, edge2D
         if (recv_from_pe(n) == prank) then
            r_count = r_count + 1
            com_edge2D%rlist(r_count) = n
         end if
      end do
   end do

   do np = 1,sPEnum
      prank = sPE(np)
      do n = 1, edge2D
         if(any(send_to_pes(:,n) == prank)) then 
            s_count = s_count + 1
            com_edge2D%slist(s_count) = n
         end if
      end do
   end do

   deallocate(send_to_pes, recv_from_pe)

end subroutine communication_edgen
!==========================================================================
subroutine mymesh(partit, mesh)
  use MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                  :: n, counter, q, k, elem, q2, eledges(4)
  integer, allocatable     :: aux(:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  !======= NODES 

  ! Owned nodes + external nodes which I need:
!!$  myDim_nod2D=count(part(:) == mype)
!!$  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1
  ! Check if the length of myList_nod2D is sufficient

!!$  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
!!$  counter=0   
!!$  do n=1, nod2D
!!$     if (part(n)==mype) then
!!$        counter=counter+1
!!$        myList_nod2D(counter)=n
!!$     end if
!!$  end do
!!$  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
!!$       com_nod2D%rlist(1:eDim_nod2D)
  ! Summary:  	     
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs

  !======= ELEMENTS
  ! 2D elements 
  ! Element belongs to PE if any of its nodes is owned by PE
  ! Element is external if it is a neighbor and does not contain owned nodes
  ! The external part is needed for FVCOM type of discretization.
  
!!$  counter=0
!!$  do n=1, elem2D
!!$     q2 = merge(3,4,elem2D_nodes(1,n) == elem2D_nodes(4,n))
!!$     do q=1,q2
!!$        if(part(elem2D_nodes(q,n))==mype) then
!!$           counter=counter+1
!!$           myList_elem2D(counter)=n
!!$           exit
!!$        end if
!!$     end do
!!$  end do
!!$  myDim_elem2D=counter
!!$  eDim_elem2D=com_elem2D%rptr(com_elem2D%rPEnum+1)-1   
  ! =======
  ! full element neighbourhood requires 
  ! a longer list     
  ! ======= 
!!$  allocate(aux(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1))    
!!$  aux=0
!!$  do n=1,com_elem2D%rptr(com_elem2D%rPEnum+1)-1         
!!$     k=com_elem2D%rlist(n)
!!$     do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$       if(com_elem2D_full%rlist(q)==k) aux(q)=1
!!$     end do
!!$  end do
!!$  ! Test: 
!!$  if(sum(aux).ne.eDim_elem2D) write(*,*) 'mymesh problem'
!!$  
!!$  counter=0
!!$  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$     if(aux(q)==0) counter=counter+1   
!!$  end do
!!$  eXDim_elem2D=counter
!!$
!!$  myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D)=&
!!$       com_elem2D%rlist(1:eDim_elem2D)
!!$  counter=myDim_elem2D+eDim_elem2D
!!$  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
!!$     if(aux(q)==0) then
!!$     counter=counter+1
!!$     myList_elem2D(counter)=com_elem2D_full%rlist(q)
!!$     end if   
!!$  end do
!!$  deallocate(aux)
              
  ! Summary: 
  ! myList_elem2D(1:myDim_elem2D) contains elements with at least single owned node.
  ! myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D) contains elements-neighbours.
  ! They are needed when interpolation is done in FV type of code.
  ! The final piece from myDim_elem2D+eDim_elem2D to 
  ! myDim_elem2D+eDim_elem2D+eXDim_elem2D is needed for MUSCL-type advection

! ======== EDGES 
  ! Owned edges (both nodes are mine)+ shared edges I do computations 
  ! at (only one node is mine; some other PE updates them
  ! simultaneously with me but I do not care ) + 
  ! external edges which I need (neither of nodes is mine, but they
  ! belong to elements in myList_elem:
  counter=0
  do n=1, edge2D
     do q=1,2 
        if (part(edges(q,n))==mype) then
         !# ??? this way an edge is added two times to the list if both if its
         ! nodes are in mype. -> no, bc of the exit statement
         ! still: why not use if (any(edges(:,n) == mype)) without
         ! a loop over the two nodes of edge n?
           counter=counter+1
           myList_edge2D(counter)=n
           exit
        end if
     end do
  end do
  myDim_edge2D=counter   ! It is the total number of my edges 
  Do n=1,myDim_elem2D
     elem=myList_elem2D(n)
     eledges=elem_edges(:,elem)
     q2 = merge(3,4,eledges(1) == eledges(4))
     ! loop over all edges of an element
     DO q=1,q2
      ! check if both nodes of the edge are not in mype and if the edge is not
      ! already collected. then add it to the local list of edges myList_edge2D
        if((part(edges(1,eledges(q))).ne.mype).and.(part(edges(2,eledges(q))).ne.mype) &
             .and. all(myList_edge2D(myDim_edge2D:counter) /= eledges(q))) then
           counter=counter+1 
           myList_edge2D(counter)=eledges(q) 
        end if
     END DO
  END DO
  eDim_edge2D=counter-myDim_edge2D
  ! Summary:  	     
  ! myList_edge2D(1:myDim_edge2D) contains owned edges +
  ! shared edges which mype updates
  ! myList_edge2D(myDim_edge2D+1:myDim_edge2D+eDim_edge2D)
  ! contains external edges which mype needs;    
end subroutine mymesh
!=================================================================
