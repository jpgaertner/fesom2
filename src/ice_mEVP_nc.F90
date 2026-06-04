module ice_mEVPdynamics_nc_interface
    interface
        subroutine mEVPdynamics_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine mEVPdynamics_div_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
   end interface
end module

module ice_mEVP_nc_interfaces
    interface
        subroutine ice_strength_mass_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine stress_tensor_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine ssh2rhs_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine stress2rhs_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine stress_tensor_div_nc(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
   end interface
end module

module mod_nc_stabilization_loc
contains
subroutine nc_stabilization_loc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_config
    use o_arrays
    use g_comm_auto
    implicit none
    type(t_ice),    intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer        :: el, eledges(3)
    real(kind=WP)  :: uu(3), vv(3), cx(3)
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice_aux, v_ice_aux, rhs_u, rhs_v
    real(kind=WP), contiguous, dimension(:), pointer  :: u_diff, v_diff
    real(kind=WP), contiguous, dimension(:), pointer  :: zeta_e
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice_aux => ice%uice_aux
    v_ice_aux => ice%vice_aux
    u_diff    => ice%nc%u_diff
    v_diff    => ice%nc%v_diff
    rhs_u     => ice%nc%rhs_u
    rhs_v     => ice%nc%rhs_v
    zeta_e    => ice%nc%zeta_e


    ! cycle 1: assemble differences
    u_diff = 0.0_WP
    v_diff = 0.0_WP

    do el = 1, myDim_elem2D
        eledges = elem_edges(:, el)

        uu = u_ice_aux(eledges)
        vv = v_ice_aux(eledges)
        u_diff(eledges(1)) = u_diff(eledges(1)) - uu(2) + uu(3)
        v_diff(eledges(1)) = v_diff(eledges(1)) - vv(2) + vv(3)
        u_diff(eledges(2)) = u_diff(eledges(2)) + uu(1) - uu(3)
        v_diff(eledges(2)) = v_diff(eledges(2)) + vv(1) - vv(3)
        u_diff(eledges(3)) = u_diff(eledges(3)) - uu(1) + uu(2)
        v_diff(eledges(3)) = v_diff(eledges(3)) - vv(1) + vv(2)
    end do

    call exchange_edge2D(u_diff, partit)
    call exchange_edge2D(v_diff, partit)

    ! cycle 2: assemble contributions to test functions
    do el = 1, myDim_elem2D
        eledges = elem_edges(:, el)

        cx = - ice%nc_stab * zeta_e(eledges)
        uu = u_diff(eledges)
        vv = v_diff(eledges)

        rhs_u(eledges(1)) = rhs_u(eledges(1)) + cx(1) * (uu(2) - uu(3))
        rhs_v(eledges(1)) = rhs_v(eledges(1)) + cx(1) * (vv(2) - vv(3))
        rhs_u(eledges(2)) = rhs_u(eledges(2)) + cx(2) * (uu(3) - uu(1))
        rhs_v(eledges(2)) = rhs_v(eledges(2)) + cx(2) * (vv(3) - vv(1))
        rhs_u(eledges(3)) = rhs_u(eledges(3)) + cx(3) * (uu(1) - uu(2))
        rhs_v(eledges(3)) = rhs_v(eledges(3)) + cx(3) * (vv(1) - vv(2))
    end do

end subroutine nc_stabilization_loc
end module

subroutine ice_strength_mass_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer        :: elem, ed, elnodes(3), eledges(3)
    real(kind=WP)  :: mmass, val3, msum, a_ice_ed, cluster_area, asum
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: ice_strength
    real(kind=WP), contiguous, dimension(:), pointer  :: zeta_e
    real(kind=WP), contiguous, dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), contiguous, dimension(:), pointer  :: inv_mass, inv_mass_a
    real(kind=WP), pointer                            :: rhoice, rhosno
    logical, contiguous, dimension(:), pointer        :: ice_exists, ice_el
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ice_strength => ice%nc%ice_strength
    zeta_e       => ice%nc%zeta_e
    a_ice        => ice%data(1)%values
    m_ice        => ice%data(2)%values
    m_snow       => ice%data(3)%values
    rhoice       => ice%thermo%rhoice
    rhosno       => ice%thermo%rhosno
    inv_mass     => ice%nc%inv_mass
    inv_mass_a   => ice%nc%inv_mass_a
    ice_exists   => ice%nc%ice_exists
    ice_el       => ice%nc%ice_el

    val3 = 1.0_WP / 3.0_WP


    zeta_e = 0.0_WP
    do elem = 1, myDim_elem2D
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1) cycle

        elnodes = elem2D_nodes(:,elem)
        eledges = elem_edges(:,elem)
        ice_el(elem) = .false.

        msum = sum(m_ice(elnodes)) * val3
        if (msum > 0.01_WP) then
            !#####
            ice_strength(elem) = ice%pstar * sum(m_ice(elnodes) * exp(-ice%c_pressure * (1.0_WP - a_ice(elnodes)))) * val3
            zeta_e(eledges) = zeta_e(eledges) + 0.5_WP * ice_strength(elem) * elem_area(elem) / ice%ice_dt
            !#####

            !##### alternative implementation
            !asum = sum(a_ice(elnodes)) * val3
            !ice_strength(elem) = ice%pstar * msum * exp(-ice%c_pressure * (1.0_WP - asum))
            !zeta_e(eledges) = zeta_e(eledges) + 0.5_WP * ice%pstar * sum(m_ice(elnodes) * exp(-ice%c_pressure*(1.0_WP-a_ice(elnodes)))) * val3 * elem_area(elem) / ice%ice_dt
            !#####
            ice_el(elem) = .true.
        end if
    end do
    call exchange_edge2D(zeta_e, partit)

    do ed = 1, myDim_edge2D
        inv_mass(ed) = 0.0_WP
        inv_mass_a(ed) = 0.0_WP
        ice_exists(ed) = .false.

        if (a_ice(edges(1,ed)) >= 0.01_WP .and. a_ice(edges(2,ed)) >= 0.01_WP) then

            mmass = 0.5_WP * sum(rhoice * m_ice(edges(:,ed)) + rhosno * m_snow(edges(:,ed)))

            !#####
            inv_mass(ed) = 1.0_WP / max(0.1, mmass)
            inv_mass_a(ed) = 0.5_WP * inv_mass(ed) * sum(a_ice(edges(:,ed)))
            !#####

            !##### alternative implementation
            !a_ice_ed = 0.5_WP * sum(a_ice(edges(:,ed)))
            !inv_mass_a(ed) = mmass / a_ice_ed
            !inv_mass_a(ed) = 1.0_WP / max(inv_mass_a(ed), 9.0_WP)

            ! the element area associated with the edge
            !if (edge_tri(2,ed)>0) then
            !    cluster_area = sum(elem_area(edge_tri(:,ed))) * val3 
            !else
            !    cluster_area = elem_area(edge_tri(1,ed)) * val3
            !end if
            !inv_mass(ed) = mmass / ((1 + mmass**2) * cluster_area)
            !#####

            ice_exists(ed) = .true.
        end if
    end do

end subroutine ice_strength_mass_nc

subroutine stress_tensor_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    implicit none 
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(inout), target  :: mesh
    !___________________________________________________________________________
    integer        :: elem, elnodes(3), eledges(3)
    real(kind=WP)  :: dx(3), dy(3), usum, vsum
    real(kind=WP)  :: val3, vale, det2, det1
    real(kind=WP)  :: x, y, meancos
    real(kind=WP)  :: eps11, eps12, eps22, eps1, eps2, pressure
    real(kind=WP)  :: r1, r2, r3, si1, si2
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice_aux, v_ice_aux
    real(kind=WP), contiguous, dimension(:), pointer  :: ice_strength
    real(kind=WP), contiguous, dimension(:), pointer  :: sigma11, sigma12, sigma22, delta
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice_aux    => ice%uice_aux
    v_ice_aux    => ice%vice_aux
    ice_strength => ice%nc%ice_strength
    sigma11      => ice%work%sigma11
    sigma12      => ice%work%sigma12
    sigma22      => ice%work%sigma22
    delta        => ice%nc%delta

    val3 = 1.0_WP / 3.0_WP
    vale = 1.0_WP / (ice%ellipse**2)
    det2 = 1.0_WP / (1.0_WP + ice%alpha_evp)
    det1 = ice%alpha_evp * det2


    do elem = 1, myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
        eledges = elem_edges(:,elem)

        dx = - 2.0_WP * gradient_sca(1:3,elem)
        dy = - 2.0_WP * gradient_sca(4:6,elem)

        ! metrics
        usum = sum(u_ice_aux(eledges))
        vsum = sum(v_ice_aux(eledges))
        meancos = metric_factor(elem)

        ! deformation rate tensor of element elem
        eps11 = sum(dx * u_ice_aux(eledges))
        eps11 = eps11 - val3 * vsum * meancos            ! metrics
        eps22 = sum(dy * v_ice_aux(eledges))
        eps12 = 0.5_WP * sum(dy * u_ice_aux(eledges) + dx * v_ice_aux(eledges))
        eps12 = eps12 + 0.5_WP * val3 * usum * meancos   ! metrics

        eps1 = eps11 + eps22
        eps2 = eps11 - eps22

        delta(elem) = eps1**2 + vale * (eps2**2 + 4.0_WP * eps12**2)
        delta(elem) = sqrt(delta(elem))
        pressure = ice_strength(elem) / (delta(elem) + ice%delta_min)

        r1 = pressure * (eps1 - delta(elem))
        r2 = pressure * eps2 * vale
        r3 = pressure * eps12 * vale
        si1 = sigma11(elem) + sigma22(elem)
        si2 = sigma11(elem) - sigma22(elem)

        si1 = det1 * si1 + det2 * r1
        si2 = det1 * si2 + det2 * r2
        sigma12(elem) = det1 * sigma12(elem) + det2 * r3
        sigma11(elem) = 0.5_WP * (si1 + si2)
        sigma22(elem) = 0.5_WP * (si1 - si2)
    end do

end subroutine stress_tensor_nc

subroutine ssh2rhs_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    implicit none 
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer        :: elem, elnodes(3), eledges(3)
    real(kind=WP)  :: dx(3), dy(3), vol
    real(kind=WP)  :: val3, aa, bb
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: gsshx, gsshy
    real(kind=WP), contiguous, dimension(:), pointer  :: elevation
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    gsshx     => ice%nc%gsshx
    gsshy     => ice%nc%gsshy
    elevation => ice%srfoce_ssh

    val3 = 1.0_WP / 3.0_WP


    ! use gsshx, gsshy for storing the contribution from elevation
    gsshx = 0.0_WP
    gsshy = 0.0_WP

    do elem = 1, myDim_elem2D
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1)  cycle

        elnodes = elem2D_nodes(:,elem)
        eledges = elem_edges(:,elem)

        vol = elem_area(elem)
        dx = gradient_sca(1:3,elem)
        dy = gradient_sca(4:6,elem)

        bb = g * val3 * vol
        aa = bb * sum(dx * elevation(elnodes))
        bb = bb * sum(dy * elevation(elnodes))

        gsshx(eledges) = gsshx(eledges) - aa
        gsshy(eledges) = gsshy(eledges) - bb
    end do
    ! integrated over triangle areas

end subroutine ssh2rhs_nc

subroutine stress2rhs_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use mod_nc_stabilization_loc
    implicit none
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(inout), target  :: mesh
    !___________________________________________________________________________
    integer        :: elem, eledges(3), k, row
    real(kind=WP)  :: val3, vol, cluster_area
    real(kind=WP)  :: x, y, meancos
    real(kind=WP)  :: dx(3), dy(3)
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: rhs_u, rhs_v
    real(kind=WP), contiguous, dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), contiguous, dimension(:), pointer  :: gsshx, gsshy, inv_mass
    logical, contiguous, dimension(:), pointer        :: ice_el
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    rhs_u    => ice%nc%rhs_u
    rhs_v    => ice%nc%rhs_v
    sigma11  => ice%work%sigma11
    sigma12  => ice%work%sigma12
    sigma22  => ice%work%sigma22
    gsshx    => ice%nc%gsshx
    gsshy    => ice%nc%gsshy
    inv_mass => ice%nc%inv_mass
    ice_el   => ice%nc%ice_el

    val3 = 1.0_WP / 3.0_WP


    rhs_u = 0.0_WP
    rhs_v = 0.0_WP

    do elem = 1, myDim_elem2D
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1)  cycle

        eledges = elem_edges(:,elem)

        vol = elem_area(elem)
        dx = - 2.0_WP * gradient_sca(1:3,elem)
        dy = - 2.0_WP * gradient_sca(4:6,elem)

        meancos = metric_factor(elem)

        if (ice_el(elem)) then
            rhs_u(eledges) = rhs_u(eledges) - vol * (sigma11(elem) * dx + sigma12(elem) * dy) &
                             - vol * sigma12(elem) * val3 * meancos     ! metrics
            rhs_v(eledges) = rhs_v(eledges) - vol * (sigma12(elem) * dx + sigma22(elem) * dy) &
                             + vol * sigma11(elem) * val3 * meancos     ! metrics
        end if
    end do

    call nc_stabilization_loc(ice, partit, mesh)

    do row = 1, myDim_edge2D
        if (myList_edge2D(row) > edge2D_in) cycle

        if (edge_tri(2,row)>0) then
            cluster_area = sum(elem_area(edge_tri(:,row))) * val3 
        else
            cluster_area = elem_area(edge_tri(1,row)) * val3
        end if

        !#####
        rhs_u(row) = ( rhs_u(row) * inv_mass(row) + gsshx(row) ) / cluster_area
        rhs_v(row) = ( rhs_v(row) * inv_mass(row) + gsshy(row) ) / cluster_area
        !#####

        !##### alternative implementation
        !rhs_u(row) = rhs_u(row) * inv_mass(row) + gsshx(row) / cluster_area
        !rhs_v(row) = rhs_v(row) * inv_mass(row) + gsshy(row) / cluster_area
        !#####
    end do

end subroutine stress2rhs_nc

subroutine mEVPdynamics_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer        :: steps, shortsteps, i
    real(kind=WP)  :: rdt, drag, det, fc
    real(kind=WP)  :: uw, vw, umod, stx, sty, rhsu, rhsv
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice, v_ice, u_ice_aux, v_ice_aux, rhs_u, rhs_v
    real(kind=WP), contiguous, dimension(:), pointer  :: inv_mass_a
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice_nod, v_ice_nod
    real(kind=WP), contiguous, dimension(:), pointer  :: u_w, v_w, stress_atmice_x, stress_atmice_y
    logical, contiguous, dimension(:), pointer        :: ice_exists
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice
    v_ice           => ice%vice
    u_ice_aux       => ice%uice_aux
    v_ice_aux       => ice%vice_aux
    rhs_u           => ice%nc%rhs_u
    rhs_v           => ice%nc%rhs_v
    u_w             => ice%srfoce_u
    v_w             => ice%srfoce_v
    stress_atmice_x => ice%stress_atmice_x
    stress_atmice_y => ice%stress_atmice_y
    inv_mass_a      => ice%nc%inv_mass_a
    u_ice_nod       => ice%nc%u_ice_nod
    v_ice_nod       => ice%nc%v_ice_nod
    ice_exists      => ice%nc%ice_exists

    steps = ice%evp_rheol_steps
    rdt = ice%ice_dt


    u_ice_aux = u_ice
    v_ice_aux = v_ice

    rhs_u = 0.0_WP
    rhs_v = 0.0_WP

    call ice_strength_mass_nc(ice, partit, mesh)
    call ssh2rhs_nc(ice, partit, mesh)

    do shortsteps = 1, steps
        call stress_tensor_nc(ice, partit, mesh)
        call stress2rhs_nc(ice, partit, mesh)

        do i = 1, myDim_edge2D
            if (myList_edge2D(i) > edge2D_in) cycle

            ! if cavity edge skip it (before checking if the second element has a cavity node, check if it exists)
            if (ulevels(edge_tri(1,i))>1) cycle
            if (edge_tri(2,i)>0) then
                if (ulevels(edge_tri(2,i))>1) cycle
            endif

            if (ice_exists(i)) then
                uw = 0.5_WP * sum(u_w(edges(:,i)))
                vw = 0.5_WP * sum(v_w(edges(:,i)))
                umod = sqrt((u_ice_aux(i) - uw)**2 + (v_ice_aux(i) - vw)**2)

                sty = inv_mass_a(i)
                drag = rdt * ice%cd_oce_ice * umod * density_0 * sty

                stx = 0.5_WP * sum(stress_atmice_x(edges(:,i))) * sty
                sty = 0.5_WP * sum(stress_atmice_y(edges(:,i))) * sty

                !rhs for water stress, air stress, and rhs_u/v (internal stress + ssh)
                rhsu = u_ice(i) + drag * uw + rdt * (stx + rhs_u(i))
                rhsv = v_ice(i) + drag * vw + rdt * (sty + rhs_v(i))

                rhsu = ice%beta_evp * u_ice_aux(i) + rhsu
                rhsv = ice%beta_evp * v_ice_aux(i) + rhsv

                !solve (Coriolis and water stress are treated implicitly)
                fc = rdt * 0.5_WP * sum(mesh%coriolis_node(edges(:,i)))
                det = (1.0_WP + ice%beta_evp + drag)**2 + fc**2
                det = minval(bc_index_nod2D(edges(:,i))) / det

                u_ice_aux(i) = det * ((1.0_WP + ice%beta_evp + drag) * rhsu + fc * rhsv)
                v_ice_aux(i) = det * ((1.0_WP + ice%beta_evp + drag) * rhsv - fc * rhsu)
            end if
        end do

        call exchange_edge2D(u_ice_aux, partit)
        call exchange_edge2D(v_ice_aux, partit)
    end do

    u_ice = u_ice_aux
    v_ice = v_ice_aux


    ! interpolate velocities from edges to nodes
    u_ice_nod = 0.0_WP
    v_ice_nod = 0.0_WP

    do i = 1, myDim_edge2D
        ! sum the adjacent edge velocities at the nodes
        u_ice_nod(edges(1,i)) = u_ice_nod(edges(1,i)) + u_ice(i)
        u_ice_nod(edges(2,i)) = u_ice_nod(edges(2,i)) + u_ice(i)
        v_ice_nod(edges(1,i)) = v_ice_nod(edges(1,i)) + v_ice(i)
        v_ice_nod(edges(2,i)) = v_ice_nod(edges(2,i)) + v_ice(i)
    end do

    do i = 1, myDim_nod2D
        ! and divide by the number of adjacent edges
        u_ice_nod(i) = u_ice_nod(i) / mesh%nn_num(i)
        v_ice_nod(i) = v_ice_nod(i) / mesh%nn_num(i)
    end do

    call exchange_nod(u_ice_nod, partit)
    call exchange_nod(v_ice_nod, partit)

end subroutine mEVPdynamics_nc

subroutine stress_tensor_div_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_comm_auto
    use mod_nc_stabilization_loc
    use ice_mEVP_nc_interfaces
    implicit none
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(inout), target  :: mesh
    !___________________________________________________________________________
    integer        :: elem, eledges(3)
    real(kind=WP)  :: val3, vale
    real(kind=WP)  :: dx(3), dy(3), usum, vsum
    real(kind=WP)  :: eps11, eps22, eps12, eps1, eps2, pressure, si11, si22, si12
    real(kind=WP)  :: meancos
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: rhs_u, rhs_v, u_ice_aux, v_ice_aux
    real(kind=WP), contiguous, dimension(:), pointer  :: ice_strength, delta
    logical, contiguous, dimension(:), pointer        :: ice_el
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    rhs_u        => ice%nc%rhs_u
    rhs_v        => ice%nc%rhs_v
    u_ice_aux    => ice%uice_aux
    v_ice_aux    => ice%vice_aux
    ice_strength => ice%nc%ice_strength
    ice_el       => ice%nc%ice_el
    delta        => ice%nc%delta

    val3 = 1.0_WP / 3.0_WP
    vale = 1.0_WP / (ice%ellipse**2)


    rhs_u = 0.0_WP
    rhs_v = 0.0_WP

    do elem = 1, myDim_elem2D
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1)  cycle

        eledges = elem_edges(:,elem)
        
        dx = - 2.0_WP * gradient_sca(1:3,elem)
        dy = - 2.0_WP * gradient_sca(4:6,elem)
        
        usum = sum(u_ice_aux(eledges))
        vsum = sum(v_ice_aux(eledges))
        meancos = metric_factor(elem)

        eps11 = sum(dx * u_ice_aux(eledges))
        eps11 = eps11 - val3 * vsum * meancos               ! metrics
        eps22 = sum(dy * v_ice_aux(eledges))
        eps12 = 0.5_WP * sum(dy * u_ice_aux(eledges) + dx * v_ice_aux(eledges))
        eps12 = eps12 + 0.5_WP * val3 * usum * meancos      ! metrics
        
        eps1 = eps11 + eps22
        eps2 = eps11 - eps22

        delta(elem) = eps1**2 + vale * (eps2**2 + 4.0_WP * eps12**2)
        delta(elem) = sqrt(delta(elem))

        pressure = ice_strength(elem) / (delta(elem) + ice%delta_min)

        si11 = 0.5_WP * pressure * (eps1 - delta(elem) + eps2 * vale)
        si22 = 0.5_WP * pressure * (eps1 - delta(elem) - eps2 * vale)
        si12 = pressure * eps12 * vale

        if (ice_el(elem)) then
            rhs_u(eledges) = rhs_u(eledges) - elem_area(elem) * (si11 * dx + si12 * dy + si12 * val3 * meancos)
            rhs_v(eledges) = rhs_v(eledges) - elem_area(elem) * (si12 * dx + si22 * dy - si11 * val3 * meancos)
        end if
    end do

    call nc_stabilization_loc(ice, partit, mesh)

end subroutine stress_tensor_div_nc

subroutine mEVPdynamics_div_nc(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_comm_auto
    use ice_mEVP_nc_interfaces
    implicit none
    type(t_ice)   , intent(inout), target  :: ice
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    !___________________________________________________________________________
    integer        :: steps, shortsteps, i
    real(kind=WP)  :: det1, det2, val3, cluster_area
    real(kind=WP)  :: rdt, drag, det, fc
    real(kind=WP)  :: uw, vw, umod, stx, sty, rhsu, rhsv
    !___________________________________________________________________________
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice, v_ice, u_ice_aux, v_ice_aux
    real(kind=WP), contiguous, dimension(:), pointer  :: rhs_u, rhs_v, R_u, R_v, gsshx, gsshy
    real(kind=WP), contiguous, dimension(:), pointer  :: inv_mass, inv_mass_a
    real(kind=WP), contiguous, dimension(:), pointer  :: u_ice_nod, v_ice_nod
    real(kind=WP), contiguous, dimension(:), pointer  :: u_w, v_w, stress_atmice_x, stress_atmice_y
    logical, contiguous, dimension(:), pointer        :: ice_exists
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice
    v_ice           => ice%vice
    u_ice_aux       => ice%uice_aux
    v_ice_aux       => ice%vice_aux
    rhs_u           => ice%nc%rhs_u
    rhs_v           => ice%nc%rhs_v
    R_u             => ice%nc%R_u
    R_v             => ice%nc%R_v
    gsshx           => ice%nc%gsshx
    gsshy           => ice%nc%gsshy
    u_w             => ice%srfoce_u
    v_w             => ice%srfoce_v
    stress_atmice_x => ice%stress_atmice_x
    stress_atmice_y => ice%stress_atmice_y
    inv_mass        => ice%nc%inv_mass
    inv_mass_a      => ice%nc%inv_mass_a
    u_ice_nod       => ice%nc%u_ice_nod
    v_ice_nod       => ice%nc%v_ice_nod
    ice_exists      => ice%nc%ice_exists

    det2 = 1.0_WP / (1.0_WP + ice%alpha_evp)
    det1 = ice%alpha_evp * det2
    val3 = 1.0_WP / 3.0_WP
    steps = ice%evp_rheol_steps
    rdt = ice%ice_dt


    call exchange_nod(bc_index_nod2D, partit)

    u_ice_aux = u_ice
    v_ice_aux = v_ice

    call ice_strength_mass_nc(ice, partit, mesh)
    call ssh2rhs_nc(ice, partit, mesh)

    do shortsteps = 1, steps
        call stress_tensor_div_nc(ice, partit, mesh)

        do i = 1, myDim_edge2D

            if (myList_edge2D(i) > edge2D_in) cycle

            ! if cavity edge skip it (before checking if the second element has a cavity node, check if it exists)
            if (ulevels(edge_tri(1,i))>1) cycle
            if (edge_tri(2,i)>0) then
                if (ulevels(edge_tri(2,i))>1) cycle
            endif

            if (ice_exists(i)) then
                ! time stepping of stress divergence
                R_u(i) = R_u(i) * det1 + det2 * rhs_u(i) * inv_mass(i)
                R_v(i) = R_v(i) * det1 + det2 * rhs_v(i) * inv_mass(i)
                if (edge_tri(2,i)>0) then
                    cluster_area = sum(elem_area(edge_tri(:,i))) * val3 
                else
                    cluster_area = elem_area(edge_tri(1,i)) * val3
                end if

                uw = 0.5_WP * sum(u_w(edges(:,i)))
                vw = 0.5_WP * sum(v_w(edges(:,i)))
                umod = sqrt((u_ice_aux(i) - uw)**2 + (v_ice_aux(i) - vw)**2)

                sty = inv_mass_a(i)
                drag = rdt * ice%cd_oce_ice * umod * density_0 * sty

                stx = 0.5_WP * sum(stress_atmice_x(edges(:,i))) * sty
                sty = 0.5_WP * sum(stress_atmice_y(edges(:,i))) * sty

                !rhs for water stress, air stress, and rhs_u/v (internal stress + ssh)
                rhsu = u_ice(i) + drag * uw + rdt * (stx + (R_u(i) + gsshx(i)) / cluster_area)
                rhsv = v_ice(i) + drag * vw + rdt * (sty + (R_v(i) + gsshy(i)) / cluster_area)

                rhsu = ice%beta_evp * u_ice_aux(i) + rhsu
                rhsv = ice%beta_evp * v_ice_aux(i) + rhsv

                !solve (Coriolis and water stress are treated implicitly)
                fc = rdt * 0.5_WP * sum(mesh%coriolis_node(edges(:,i)))
                det = (1.0_WP + ice%beta_evp + drag)**2 + fc**2
                det = minval(bc_index_nod2D(edges(:,i))) / det

                u_ice_aux(i) = det * ((1.0_WP + ice%beta_evp + drag) * rhsu + fc * rhsv)
                v_ice_aux(i) = det * ((1.0_WP + ice%beta_evp + drag) * rhsv - fc * rhsu)
            end if
        end do

        call exchange_edge2D(u_ice_aux, partit)
        call exchange_edge2D(v_ice_aux, partit)
    end do

    u_ice = u_ice_aux
    v_ice = v_ice_aux


    ! interpolate velocities from edges to nodes
    u_ice_nod = 0.0_WP
    v_ice_nod = 0.0_WP

    do i = 1, myDim_edge2D
        ! sum the adjacent edge velocities at the nodes
        u_ice_nod(edges(1,i)) = u_ice_nod(edges(1,i)) + u_ice(i)
        u_ice_nod(edges(2,i)) = u_ice_nod(edges(2,i)) + u_ice(i)
        v_ice_nod(edges(1,i)) = v_ice_nod(edges(1,i)) + v_ice(i)
        v_ice_nod(edges(2,i)) = v_ice_nod(edges(2,i)) + v_ice(i)
    end do

    do i = 1, myDim_nod2D
        ! and divide by the number of adjacent edges
        u_ice_nod(i) = u_ice_nod(i) / mesh%nn_num(i)
        v_ice_nod(i) = v_ice_nod(i) / mesh%nn_num(i)
    end do

    call exchange_nod(u_ice_nod, partit)
    call exchange_nod(v_ice_nod, partit)

end subroutine mEVPdynamics_div_nc
