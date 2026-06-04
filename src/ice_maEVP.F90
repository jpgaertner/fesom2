module ice_maEVP_interfaces
    interface
        subroutine ssh2rhs(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine ssh2rhs

        subroutine stress_tensor_a(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine stress_tensor_a

        subroutine stress2rhs_m(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine stress2rhs_m

        subroutine find_alpha_field_a(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine find_alpha_field_a

        subroutine find_beta_field_a(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine find_beta_field_a
   end interface
end module ice_maEVP_interfaces

module ice_maEVPdynamics_interface
    interface
        subroutine EVPdynamics_a(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_mesh)  , intent(in)   , target :: mesh
        type(t_partit), intent(inout), target :: partit
        type(t_ice)   , intent(inout), target :: ice
        end subroutine EVPdynamics_a

        subroutine EVPdynamics_m(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_mesh)  , intent(in)   , target :: mesh
        type(t_partit), intent(inout), target :: partit
        type(t_ice)   , intent(inout), target :: ice
        end subroutine EVPdynamics_m
   end interface
end module ice_maEVPdynamics_interface

module mod_nc_stabilization
contains
subroutine nc_stabilization(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_config
    use o_arrays
    use g_comm_auto
    implicit none
    type(t_ice),    intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: el, eledges(3)
    real(kind=WP)                         :: uu(3), vv(3), cx(3)
    real(kind=WP), allocatable            :: udif(:), vdif(:)
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux, u_rhs_ice, v_rhs_ice, zeta_e
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice_aux => ice%uice_aux(:)
    v_ice_aux => ice%vice_aux(:)
    u_rhs_ice => ice%uice_rhs(:)
    v_rhs_ice => ice%vice_rhs(:)
    zeta_e    => ice%zeta_e(:)

    ! cycle 1: assemble differences
    allocate(udif(myDim_edge2D+eDim_edge2D), vdif(myDim_edge2D+eDim_edge2D))
    udif = 0.0_WP
    vdif = 0.0_WP

    do el = 1, myDim_elem2D
        eledges = elem_edges(:, el)

        uu = u_ice_aux(eledges)
        vv = v_ice_aux(eledges)
        udif(eledges(1))=udif(eledges(1))-uu(2)+uu(3)
        vdif(eledges(1))=vdif(eledges(1))-vv(2)+vv(3)
        udif(eledges(2))=udif(eledges(2))+uu(1)-uu(3)
        vdif(eledges(2))=vdif(eledges(2))+vv(1)-vv(3)
        udif(eledges(3))=udif(eledges(3))-uu(1)+uu(2)
        vdif(eledges(3))=vdif(eledges(3))-vv(1)+vv(2)
    end do

    call exchange_edge2D(udif, partit)
    call exchange_edge2D(vdif, partit)

    ! cycle 2: assemble contributions to test functions
    do el = 1, myDim_elem2D
        eledges = elem_edges(:, el)

        cx = -ice%nc_stab * zeta_e(eledges)
        uu = udif(eledges)
        vv = vdif(eledges)

        u_rhs_ice(eledges(1))=u_rhs_ice(eledges(1))+cx(1)*(uu(2)-uu(3))
        v_rhs_ice(eledges(1))=v_rhs_ice(eledges(1))+cx(1)*(vv(2)-vv(3))
        u_rhs_ice(eledges(2))=u_rhs_ice(eledges(2))+cx(2)*(uu(3)-uu(1))
        v_rhs_ice(eledges(2))=v_rhs_ice(eledges(2))+cx(2)*(vv(3)-vv(1))
        u_rhs_ice(eledges(3))=u_rhs_ice(eledges(3))+cx(3)*(uu(1)-uu(2))
        v_rhs_ice(eledges(3))=v_rhs_ice(eledges(3))+cx(3)*(vv(1)-vv(2))
    end do

    deallocate(udif, vdif)

end subroutine
end module mod_nc_stabilization
!
!
!_______________________________________________________________________________
! New evp implementation following Bouillion et al. 2013
! and Kimmritz et al. 2015 (mEVP) and Kimmritz et al. 2016 (aEVP)
! Internal stress tensor
! New implementation following Boullion et al, Ocean Modelling 2013.
! SD, 30.07.2014
subroutine stress_tensor_m(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use mod_mesh
    use g_config
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer         :: elem, elnodes(3)
    real(kind=WP)   :: dx(3), dy(3), msum, asum
    real(kind=WP)   :: eps1, eps2, pressure, delta
    real(kind=WP)   :: val3, meancos, usum, vsum, vale
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    eps11        => ice%work%eps11(:)
    eps12        => ice%work%eps12(:)
    eps22        => ice%work%eps22(:)
    sigma11      => ice%work%sigma11(:)
    sigma12      => ice%work%sigma12(:)
    sigma22      => ice%work%sigma22(:)
    u_ice_aux    => ice%uice_aux(:)
    v_ice_aux    => ice%vice_aux(:)

    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
    vale=1.0_WP/(ice%ellipse**2)
    det2=1.0_WP/(1.0_WP+ice%alpha_evp)
    det1=ice%alpha_evp*det2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, dx, dy, msum, asum, eps1, eps2, pressure, delta, meancos, usum, vsum, r1, r2, r3, si1, si2)
    do elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1) cycle

        msum=sum(m_ice(elnodes))*val3
        if(msum<=0.01_WP) cycle !DS
        asum=sum(a_ice(elnodes))*val3

        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        ! METRICS:
            vsum=sum(v_ice_aux(elnodes))
            usum=sum(u_ice_aux(elnodes))
            meancos=metric_factor(elem)
        !
        ! ====== Deformation rate tensor on element elem:
        eps11(elem)=sum(dx*u_ice_aux(elnodes))
        eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
        eps22(elem)=sum(dy*v_ice_aux(elnodes))
        eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
        eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics

        ! ======= Switch to eps1,eps2
        eps1=eps11(elem)+eps22(elem)
        eps2=eps11(elem)-eps22(elem)

        ! ====== moduli:
        delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
        delta=sqrt(delta)

#if defined (__icepack)
        pressure = sum(strength(elnodes))*val3/max(delta,ice%delta_min)
#else
        pressure=ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))/max(delta,ice%delta_min)
#endif
        r1=pressure*(eps1-max(delta,ice%delta_min))
        r2=pressure*eps2*vale
        r3=pressure*eps12(elem)*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_WP*(si1+si2)
        sigma22(elem)=0.5_WP*(si1-si2)

#if defined (__icepack)
        rdg_conv_elem(elem)  = -min((eps11(elem)+eps22(elem)),0.0_WP)
        rdg_shear_elem(elem) = 0.5_WP*(delta - abs(eps11(elem)+eps22(elem)))
#endif
    end do
!$OMP END PARALLEL DO
    ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of
    ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
    ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
    ! and similarly for si2 and sigma12
end subroutine stress_tensor_m

!
!
!_______________________________________________________________________________
! Compute the contribution from the elevation to the rhs
! S.D. 30.07.2014
subroutine ssh2rhs(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_param
    use mod_mesh
    use g_config
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: row, elem, elnodes(3), n
    real(kind=WP)            :: dx(3), dy(3), vol
    real(kind=WP)            :: val3, meancos, aa, bb, p_ice(3)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: rhs_a, rhs_m
    real(kind=WP), dimension(:), pointer  :: elevation
    real(kind=WP)              , pointer  :: rhoice, rhosno, inv_rhowat
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    m_ice        => ice%data(2)%values(:)
    m_snow       => ice%data(3)%values(:)
    rhs_a        => ice%data(1)%values_rhs(:)
    rhs_m        => ice%data(2)%values_rhs(:)
    elevation    => ice%srfoce_ssh
    rhoice       => ice%thermo%rhoice
    rhosno       => ice%thermo%rhosno
    inv_rhowat   => ice%thermo%inv_rhowat

    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(row, elem, elnodes, n, dx, dy, vol, meancos, aa, bb, p_ice)
!$OMP DO
    ! use rhs_m and rhs_a for storing the contribution from elevation:
    do row=1, myDim_nod2d
        rhs_a(row)=0.0_WP
        rhs_m(row)=0.0_WP
    end do
!$OMP END DO
    !_____________________________________________________________________________
    ! use floating sea ice for zlevel and zstar
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') then
!$OMP DO
        do elem=1,myDim_elem2d
            elnodes=elem2D_nodes(:,elem)
            !_______________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(elem) > 1) cycle

            !_______________________________________________________________________
            vol=elem_area(elem)
            dx=gradient_sca(1:3,elem)
            dy=gradient_sca(4:6,elem)

            !_______________________________________________________________________
            ! add pressure gradient from sea ice --> in case of floating sea ice
            p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
            do n=1,3
                p_ice(n)=min(p_ice(n),max_ice_loading)
            end do

            !_______________________________________________________________________
            bb=g*val3*vol
            aa=bb*sum(dx*(elevation(elnodes)+p_ice))
            bb=bb*sum(dy*(elevation(elnodes)+p_ice))
            do n=1,3
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(elnodes(n)))
#else
!$OMP ORDERED
#endif
               rhs_a(elnodes(n))=rhs_a(elnodes(n))-aa
               rhs_m(elnodes(n))=rhs_m(elnodes(n))-bb
#if defined(_OPENMP) && !defined(__openmp_reproducible)
               call omp_unset_lock(partit%plock(elnodes(n)))
#else
!$OMP END ORDERED
#endif
            end do
        end do
!$OMP END DO
    else
!$OMP DO
        do elem=1,myDim_elem2d
            elnodes=elem2D_nodes(:,elem)
            !_______________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(elem) > 1) cycle

            vol=elem_area(elem)
            dx=gradient_sca(1:3,elem)
            dy=gradient_sca(4:6,elem)
            bb=g*val3*vol
            aa=bb*sum(dx*elevation(elnodes))
            bb=bb*sum(dy*elevation(elnodes))
            do n=1,3
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(elnodes(n)))
#else
!$OMP ORDERED
#endif
               rhs_a(elnodes(n))=rhs_a(elnodes(n))-aa
               rhs_m(elnodes(n))=rhs_m(elnodes(n))-bb
#if defined(_OPENMP) && !defined(__openmp_reproducible)
               call omp_unset_lock(partit%plock(elnodes(n)))
#else
!$OMP END ORDERED
#endif
            end do
        end do
!$OMP END DO
    end if
!$OMP END PARALLEL
end subroutine ssh2rhs
!
!
!_______________________________________________________________________________
! add internal stress to the rhs
! SD, 30.07.2014
subroutine stress2rhs_m(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_param
    use mod_mesh
    use g_config
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: k, row, elem, elnodes(3)
    real(kind=WP)            :: dx(3), dy(3), vol
    real(kind=WP)            :: val3, mf, aa, bb
    real(kind=WP)            :: mass, cluster_area, elevation_elem(3)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
    real(kind=WP)              , pointer  :: rhoice, rhosno
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    m_snow       => ice%data(3)%values(:)
    sigma11      => ice%work%sigma11(:)
    sigma12      => ice%work%sigma12(:)
    sigma22      => ice%work%sigma22(:)
    u_rhs_ice    => ice%uice_rhs(:)
    v_rhs_ice    => ice%vice_rhs(:)
    rhs_a        => ice%data(1)%values_rhs(:)
    rhs_m        => ice%data(2)%values_rhs(:)
    rhoice       => ice%thermo%rhoice
    rhosno       => ice%thermo%rhosno

    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
!$OMP PARALLEL DO
    do row=1, myDim_nod2d
        u_rhs_ice(row)=0.0_WP
        v_rhs_ice(row)=0.0_WP
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem, elnodes, k, row, dx, dy, vol, mf, aa, bb, mass, cluster_area, elevation_elem)
!$OMP DO
    do elem=1,myDim_elem2d
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1) cycle

        if(sum(a_ice(elnodes)) < 0.01_WP) cycle !DS

        vol=elem_area(elem)
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        mf=metric_factor(elem)                               !metrics

        do k=1,3
            row=elnodes(k)
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(row))
#else
!$OMP ORDERED
#endif
            u_rhs_ice(row)=u_rhs_ice(row) - vol* &
                (sigma11(elem)*dx(k)+sigma12(elem)*dy(k))    &
        -vol*sigma12(elem)*val3*mf                         !metrics
            v_rhs_ice(row)=v_rhs_ice(row) - vol* &
                (sigma12(elem)*dx(k)+sigma22(elem)*dy(k))    &
        +vol*sigma11(elem)*val3*mf                         ! metrics
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(row))
#else
!$OMP END ORDERED
#endif
        end do
    end do
!$OMP END DO
!$OMP DO
    do row=1, myDim_nod2d
        !_______________________________________________________________________
        ! if cavity node skip it
        if ( ulevels_nod2d(row)>1 ) cycle

        mass=(m_ice(row)*rhoice+m_snow(row)*rhosno)
        mass=mass/(1.0_WP+mass*mass)
        u_rhs_ice(row)=(u_rhs_ice(row)*mass + rhs_a(row))/area(1,row)
        v_rhs_ice(row)=(v_rhs_ice(row)*mass + rhs_m(row))/area(1,row)
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine stress2rhs_m
!
!
!_______________________________________________________________________________
! assemble rhs and solve for ice velocity
! New implementation based on Bouillion et al. Ocean Modelling 2013
! SD 30.07.14
subroutine EVPdynamics_m(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_config
    use o_arrays
    use g_comm_auto
    use mod_nc_stabilization
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
    use icedrv_main,   only: icepack_to_fesom
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer          :: steps, shortstep, i, ed,n
    real(kind=WP)    :: rdt, drag, det
    real(kind=WP)    :: umod, rhsu, rhsv
    logical          :: ice_el(partit%myDim_elem2D)
    !NR for stress_tensor_m
    integer         :: el, elnodes(3), eledges(3)
    real(kind=WP)   :: msum, asum
    real(kind=WP)   :: eps1, eps2, pressure, pressure_fac(partit%myDim_elem2D)!, delta
    real(kind=WP)   :: val3, meancos, vale
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2
    !NR for stress2rhs_m
    integer        :: k, row
    real(kind=WP)  :: vol
    real(kind=WP)  :: mf,aa, bb,p_ice(3)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
    real(kind=WP), dimension(:), pointer  :: u_w, v_w
    real(kind=WP), dimension(:), pointer  :: elevation
    real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux
    real(kind=WP), dimension(:), pointer  :: ice_strength
#if defined (__icepack)
    real(kind=WP), dimension(:), pointer  :: a_ice_old, m_ice_old, m_snow_old
#endif
    real(kind=WP)              , pointer  :: rhoice, rhosno, inv_rhowat
    !___________________________________________________________________________
    integer, pointer                                :: ice_vplace
    real(kind=WP), dimension(:), pointer            :: u_ice_nod, v_ice_nod
    real(kind=WP), dimension(3,partit%myDim_elem2D) :: dx, dy
    integer, dimension(3,partit%myDim_elem2D)       :: placement
    integer                                         :: myDim, eDim
    logical, allocatable                            :: ice_exists(:)
    real(kind=WP), allocatable                      :: inv_thickness(:), mass(:)
    real(kind=WP)                                   :: a_ice_ed, area_ed, uw, vw, stx, sty, cor
    real(kind=WP), dimension(:), pointer            :: zeta_e, delta
    !___________________________________________________________________________
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    ice_vplace      => ice%ice_vplace
    u_ice           => ice%uice(:)
    v_ice           => ice%vice(:)
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    m_snow          => ice%data(3)%values(:)
    eps11           => ice%work%eps11(:)
    eps12           => ice%work%eps12(:)
    eps22           => ice%work%eps22(:)
    sigma11         => ice%work%sigma11(:)
    sigma12         => ice%work%sigma12(:)
    sigma22         => ice%work%sigma22(:)
    u_rhs_ice       => ice%uice_rhs(:)
    v_rhs_ice       => ice%vice_rhs(:)
    rhs_a           => ice%data(1)%values_rhs(:)
    rhs_m           => ice%data(2)%values_rhs(:)
    u_w             => ice%srfoce_u(:)
    v_w             => ice%srfoce_v(:)
    elevation       => ice%srfoce_ssh(:)
    stress_atmice_x => ice%stress_atmice_x(:)
    stress_atmice_y => ice%stress_atmice_y(:)
    u_ice_aux       => ice%uice_aux(:)
    v_ice_aux       => ice%vice_aux(:)
    ice_strength    => ice%work%ice_strength(:)
#if defined (__icepack)
    a_ice_old       => ice%data(1)%values_old(:)
    m_ice_old       => ice%data(2)%values_old(:)
    m_snow_old      => ice%data(3)%values_old(:)
#endif
    rhoice          => ice%thermo%rhoice
    rhosno          => ice%thermo%rhosno
    inv_rhowat      => ice%thermo%inv_rhowat
    u_ice_nod       => ice%nc%u_ice_nod(:)
    v_ice_nod       => ice%nc%v_ice_nod(:)
    zeta_e          => ice%zeta_e(:)
    delta           => ice%nc%delta(:)
    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
    vale=1.0_WP/(ice%ellipse**2)
    det2=1.0_WP/(1.0_WP+ice%alpha_evp)
    det1=ice%alpha_evp*det2
    rdt=ice%ice_dt
    steps=ice%evp_rheol_steps

    if (ice_vplace == 0) then
        do el = 1, myDim_elem2D
            dx(:,el) = gradient_sca(1:3,el)
            dy(:,el) = gradient_sca(4:6,el)
            placement(:,el) = elem2D_nodes(:,el)
        end do
        myDim = myDim_nod2D
        eDim = eDim_nod2D
    else if (ice_vplace == 1) then
        do el = 1, myDim_elem2D
            dx(:,el) = - 2.0_WP * gradient_sca(1:3,el)
            dy(:,el) = - 2.0_WP * gradient_sca(4:6,el)
            placement(:,el) = elem_edges(:,el)
        end do
        myDim = myDim_edge2D
        eDim = eDim_edge2D
    end if

    allocate(ice_exists(myDim))
    allocate(inv_thickness(myDim))
    allocate(mass(myDim))


    !___________________________________________________________________________
    u_ice_aux=u_ice    ! Initialize solver variables
    v_ice_aux=v_ice

    !___________________________________________________________________________
    ! Populate the nodal ice_strength diagnostic with the canonical Hibler
    ! (1979) P = P*·h·exp(-C(1-A)) from per-node m_ice and a_ice. Exposed as
    ! the 'strength_ice' output stream. mEVP evaluates the per-element
    ! strength inline in its iteration, so this nodal field is purely a
    ! diagnostic.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n)
    do n = 1, myDim_nod2D
        ice_strength(n) = 0.0_WP
        if (ulevels_nod2D(n) > 1) cycle
        if (m_ice(n) <= 0._WP .or. a_ice(n) <= 0._WP) cycle
        ice_strength(n) = ice%pstar*m_ice(n)*exp(-ice%c_pressure*(1.0_WP-a_ice(n)))
    end do
!$OMP END PARALLEL DO

#if defined (__icepack)
    a_ice_old(:)  = a_ice(:)
    m_ice_old(:)  = a_ice(:)
    m_snow_old(:) = m_snow(:)

    call icepack_to_fesom (nx_in=(myDim_nod2D+eDim_nod2D), &
                            aice_out=a_ice,                 &
                            vice_out=m_ice,                 &
                            vsno_out=m_snow)
#endif

    !NR inlined, to have all initialization in one place.
    !  call ssh2rhs
    ! use rhs_m and rhs_a for storing the contribution from elevation:
!$OMP PARALLEL DO
    do row=1, myDim
        rhs_a(row)=0.0_WP
        rhs_m(row)=0.0_WP
    end do
!$OMP END PARALLEL DO
    !_____________________________________________________________________________
    ! use floating sea ice for zlevel and zstar
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el, elnodes, vol, dx, dy, p_ice, n, bb, aa)
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') then
!$OMP DO
        do el=1,myDim_elem2d
            elnodes=elem2D_nodes(:,el)

            !_______________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(el) > 1) cycle

            !_______________________________________________________________________
            vol=elem_area(el)

            !_______________________________________________________________________
            ! add pressure gradient from sea ice --> in case of floating sea ice
            p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
            do n=1,3
                p_ice(n)=min(p_ice(n),max_ice_loading)
            end do

            !_______________________________________________________________________
            bb=g*val3*vol
            aa=bb*sum(gradient_sca(1:3,el)*(elevation(elnodes)+p_ice))
            bb=bb*sum(gradient_sca(4:6,el)*(elevation(elnodes)+p_ice))
            do n=1, 3
#if defined(_OPENMP) && !defined(__openmp_reproducible)
               call omp_set_lock  (partit%plock(elnodes(n)))
#else
!$OMP ORDERED
#endif
               rhs_a(placement(n,el))=rhs_a(placement(n,el)) - aa
               rhs_m(placement(n,el))=rhs_m(placement(n,el)) - bb
#if defined(_OPENMP) && !defined(__openmp_reproducible)
               call omp_unset_lock(partit%plock(elnodes(n)))
#else
!$OMP END ORDERED
#endif
            end do
        end do
!$OMP END DO
    !_____________________________________________________________________________
    ! use levitating sea ice for linfs, zlevel and zstar
    else
!$OMP DO
        do el=1,myDim_elem2d
            elnodes=elem2D_nodes(:,el)
            !_______________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(el) > 1)  cycle

            vol=elem_area(el)
            bb=g*val3*vol
            aa=bb*sum(gradient_sca(1:3,el)*elevation(elnodes))
            bb=bb*sum(gradient_sca(4:6,el)*elevation(elnodes))
            do n=1, 3
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(elnodes(n)))
#else
!$OMP ORDERED
#endif
               rhs_a(placement(n,el))=rhs_a(placement(n,el)) - aa
               rhs_m(placement(n,el))=rhs_m(placement(n,el)) - bb
#if defined(_OPENMP) && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(elnodes(n)))
#else
!$OMP END ORDERED
#endif
            end do
        end do
!$OMP END DO
    end if
!$OMP END PARALLEL
    !___________________________________________________________________________
    ! precompute thickness (the inverse is needed) and mass (scaled by area)
!$OMP PARALLEL DO
    if (ice_vplace == 0) then
        do i=1,myDim_nod2D
            inv_thickness(i) = 0._WP
            mass(i) = 0._WP
            ice_exists(i) = .false.
            !_______________________________________________________________________
            ! if cavity node skip it
            if ( ulevels_nod2d(i)>1 ) cycle

            if (a_ice(i) >= 0.01_WP) then
                inv_thickness(i) = (rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
                inv_thickness(i) = 1.0_WP/max(inv_thickness(i), 9.0_WP)  ! Limit the mass

                mass(i) = (m_ice(i)*rhoice+m_snow(i)*rhosno)
                mass(i) = mass(i)/((1.0_WP+mass(i)*mass(i))*area(1,i))

                ! scale rhs_a, rhs_m, too.
                rhs_a(i) = rhs_a(i)/area(1,i)
                rhs_m(i) = rhs_m(i)/area(1,i)

                ice_exists(i) = .true.
            endif
        enddo
    else if (ice_vplace == 1) then
        do i=1,myDim_edge2D
            inv_thickness(i) = 0.0_WP
            mass(i) = 0.0_WP
            ice_exists(i) = .false.

            ! if cavity edge skip it (before checking if the second element has a cavity node, check if it exists)
            if (ulevels(edge_tri(1,i))>1) cycle
            if (edge_tri(2,i)>0) then
                if (ulevels(edge_tri(2,i))>1) cycle
            endif

            if (a_ice(edges(1,i)) >= 0.01_WP .and. a_ice(edges(2,i)) >= 0.01_WP) then
                a_ice_ed = 0.5_WP * sum(a_ice(edges(:,i))) ! average the two adjacent nodes of the edge
                inv_thickness(i) = 0.5_WP * sum(rhoice*m_ice(edges(:,i))+rhosno*m_snow(edges(:,i)))/a_ice_ed
                inv_thickness(i) = 1.0_WP/max(inv_thickness(i), 9.0_WP)  ! Limit the mass

                ! the area associated with the edge
                if (edge_tri(2,i)>0) then
                    area_ed = sum(elem_area(edge_tri(:,i))) * val3 
                else
                    area_ed = elem_area(edge_tri(1,i)) * val3
                end if
                
                mass(i) = 0.5_WP * sum(m_ice(edges(:,i))*rhoice+m_snow(edges(:,i))*rhosno)
                mass(i) = mass(i)/((1.0_WP+mass(i)*mass(i))*area_ed)

                ! scale rhs_a, rhs_m, too.
                rhs_a(i) = rhs_a(i)/area_ed
                rhs_m(i) = rhs_m(i)/area_ed

                ice_exists(i) = .true.
            endif
        enddo
    end if

!$OMP END PARALLEL DO
    !___________________________________________________________________________
    ! precompute pressure factor
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el, elnodes, msum, asum)
            
    zeta_e = 0.0_WP
    do el=1,myDim_elem2D
        elnodes=elem2D_nodes(:,el)
        eledges=elem_edges(:,el)
        pressure_fac(el) = 0._WP
        ice_el(el) = .false.

        !_______________________________________________________________________
        ! if element has any cavity node skip it
        if (ulevels(el) > 1) cycle

        msum=sum(m_ice(elnodes))*val3
        if(msum > 0.01) then
            ice_el(el) = .true.
            asum=sum(a_ice(elnodes))*val3
            pressure_fac(el) = det2*ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
            
            zeta_e(eledges) = zeta_e(eledges) + 0.5_WP * ice%pstar * sum(m_ice(elnodes) * exp(-ice%c_pressure*(1.0_WP-a_ice(elnodes)))) * val3 * elem_area(el) / ice%ice_dt 
        endif

    end do
    call exchange_edge2D(zeta_e, partit)

!$OMP END PARALLEL DO
!$OMP PARALLEL DO
    do row=1, myDim
        u_rhs_ice(row)=0.0_WP
        v_rhs_ice(row)=0.0_WP
    end do
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    ! Ice EVPdynamics Iteration main loop:
#if defined (__icepack)
    rdg_conv_elem(:)  = 0.0_WP
    rdg_shear_elem(:) = 0.0_WP
#endif
    do shortstep=1, steps
        !NR inlining, to make it easier to have local arrays and fuse loops
        !NR    call stress_tensor_m
        ! Internal stress tensor
        ! New implementation following Boullion et al, Ocean Modelling 2013.
        ! SD, 30.07.2014
        !_______________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el, i, ed, row, elnodes, dx, dy, meancos, eps1, eps2, delta, pressure, umod, drag, rhsu, rhsv, det, n)
!$OMP DO
        do el=1,myDim_elem2D
            if (ulevels(el)>1) cycle

            !___________________________________________________________________
            if(ice_el(el)) then

                ! METRICS:
                meancos = val3*metric_factor(el)
                !
                ! ====== Deformation rate tensor on element elem:
                eps11(el) = sum(dx(:,el)*u_ice_aux(placement(:,el))) - sum(v_ice_aux(placement(:,el)))*meancos                !metrics
                eps22(el) = sum(dy(:,el)*v_ice_aux(placement(:,el)))
                eps12(el) = 0.5_WP*(sum(dy(:,el)*u_ice_aux(placement(:,el)) + dx(:,el)*v_ice_aux(placement(:,el))) &
                                +sum(u_ice_aux(placement(:,el)))*meancos )          !metrics

                ! ======= Switch to eps1,eps2
                eps1 = eps11(el) + eps22(el)
                eps2 = eps11(el) - eps22(el)

                ! ====== moduli:
                delta(el) = sqrt(eps1**2+vale*(eps2**2+4.0_WP*eps12(el)**2))

                pressure = pressure_fac(el)/(delta(el)+ice%delta_min)

                !        si1 = det1*(sigma11(el)+sigma22(el)) + pressure*(eps1-delta)
                !        si2 = det1*(sigma11(el)-sigma22(el)) + pressure*eps2*vale
                !        sigma11(el) = 0.5_WP*(si1+si2)
                !        sigma22(el) = 0.5_WP*(si1-si2)
                !NR directly insert si1, si2 cancels some operations and should increase accuracy
                sigma12(el) = det1*sigma12(el) +       pressure*eps12(el)*vale
                sigma11(el) = det1*sigma11(el) + 0.5_WP*pressure*(eps1 - delta(el) + eps2*vale)
                sigma22(el) = det1*sigma22(el) + 0.5_WP*pressure*(eps1 - delta(el) - eps2*vale)

#if defined (__icepack)
                rdg_conv_elem(el)  = -min((eps11(el)+eps22(el)),0.0_WP)
                rdg_shear_elem(el) = 0.5_WP*(delta(el) - abs(eps11(el)+eps22(el)))
#endif

                !  end do   ! fuse loops
                ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of
                ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
                ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
                ! and similarly for si2 and sigma12

                !NR inlining  call stress2rhs_m
                ! add internal stress to the rhs
                ! SD, 30.07.2014
                !-----------------------------------------------------------------
                if (placement(1,el) <= myDim) then
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_set_lock  (partit%plock(placement(1,el)))
#else
!$OMP ORDERED
#endif
                    u_rhs_ice(placement(1,el)) = u_rhs_ice(placement(1,el)) - elem_area(el) * &
                            (sigma11(el) * dx(1,el) + sigma12(el) * dy(1,el) + sigma12(el) * meancos)     !metrics
                    v_rhs_ice(placement(1,el)) = v_rhs_ice(placement(1,el)) - elem_area(el) * &
                            (sigma12(el) * dx(1,el) + sigma22(el) * dy(1,el) - sigma11(el) * meancos)     !metrics
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(placement(1,el)))
#else
!$OMP END ORDERED
#endif
                end if

                if (placement(2,el) <= myDim) then
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_set_lock  (partit%plock(placement(2,el)))
#else
!$OMP ORDERED
#endif
                    u_rhs_ice(placement(2,el)) = u_rhs_ice(placement(2,el)) - elem_area(el) * &
                            (sigma11(el) * dx(2,el) + sigma12(el) * dy(2,el) + sigma12(el) * meancos)     !metrics
                    v_rhs_ice(placement(2,el)) = v_rhs_ice(placement(2,el)) - elem_area(el) * &
                            (sigma12(el) * dx(2,el) + sigma22(el) * dy(2,el) - sigma11(el) * meancos)     !metrics
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(placement(2,el)))
#else
!$OMP END ORDERED
#endif
                end if

                if (placement(3,el) <= myDim) then
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_set_lock  (partit%plock(placement(3,el)))
#else
!$OMP ORDERED
#endif
                    u_rhs_ice(placement(3,el)) = u_rhs_ice(placement(3,el)) - elem_area(el) * &
                            (sigma11(el) * dx(3,el) + sigma12(el) * dy(3,el) + sigma12(el) * meancos)      !metrics
                    v_rhs_ice(placement(3,el)) = v_rhs_ice(placement(3,el)) - elem_area(el) * &
                            (sigma12(el) * dx(3,el) + sigma22(el) * dy(3,el) - sigma11(el) * meancos)      !metrics
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                   call omp_unset_lock(partit%plock(placement(3,el)))
#else
!$OMP END ORDERED
#endif
                end if
            end if
        end do ! --> do el=1,myDim_elem2D
!$OMP END DO
        if (ice_vplace == 1) then
            call nc_stabilization(ice, partit, mesh)
        end if
!$OMP DO
        if (ice_vplace == 0) then
            do i=1, myDim_nod2d
                !___________________________________________________________________
                if (ulevels_nod2D(i)>1) cycle

                !___________________________________________________________________
                if (ice_exists(i)) then                   ! Skip if ice is absent

                    u_rhs_ice(i) = u_rhs_ice(i)*mass(i) + rhs_a(i)
                    v_rhs_ice(i) = v_rhs_ice(i)*mass(i) + rhs_m(i)
                    ! end do   !NR fuse loops
                    !============= stress2rhs_m ends ======================
                    !    do i=1,myDim_nod2D
                    umod = sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
                    drag = rdt*ice%cd_oce_ice*umod*density_0*inv_thickness(i)

                    !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
                    rhsu = u_ice(i)+drag*u_w(i)+rdt*(inv_thickness(i)*stress_atmice_x(i)+u_rhs_ice(i)) + ice%beta_evp*u_ice_aux(i)
                    rhsv = v_ice(i)+drag*v_w(i)+rdt*(inv_thickness(i)*stress_atmice_y(i)+v_rhs_ice(i)) + ice%beta_evp*v_ice_aux(i)

                    !solve (Coriolis and water stress are treated implicitly)
                    det = bc_index_nod2D(i) / ((1.0_WP+ice%beta_evp+drag)**2 + (rdt*mesh%coriolis_node(i))**2)

                    u_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsu +rdt*mesh%coriolis_node(i)*rhsv)
                    v_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsv -rdt*mesh%coriolis_node(i)*rhsu)
                end if
            end do ! --> do i=1, myDim_nod2d
        else if (ice_vplace == 1) then
            do i=1, myDim_edge2D

                ! if cavity edge skip it (before checking if the second element has a cavity node, check if it exists)
                if (ulevels(edge_tri(1,i))>1) cycle
                if (edge_tri(2,i)>0) then
                    if (ulevels(edge_tri(2,i))>1) cycle
                endif
                
                if (ice_exists(i)) then

                    u_rhs_ice(i) = u_rhs_ice(i)*mass(i) + rhs_a(i)
                    v_rhs_ice(i) = v_rhs_ice(i)*mass(i) + rhs_m(i)

                    ! interpolate ocean surface velocity and wind stress from nodes to edges
                    uw = 0.5_WP * sum(u_w(edges(:,i)))
                    vw = 0.5_WP * sum(v_w(edges(:,i)))
                    stx = 0.5_WP * sum(stress_atmice_x(edges(:,i)))
                    sty = 0.5_WP * sum(stress_atmice_y(edges(:,i)))
                    cor = 0.5_WP * sum(mesh%coriolis_node(edges(:,i)))

                    umod = sqrt((u_ice_aux(i)-uw)**2+(v_ice_aux(i)-vw)**2)
                    drag = rdt*ice%cd_oce_ice*umod*density_0*inv_thickness(i)

                    !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
                    rhsu = u_ice(i)+drag*uw+rdt*(inv_thickness(i)*stx+u_rhs_ice(i)) + ice%beta_evp*u_ice_aux(i)
                    rhsv = v_ice(i)+drag*vw+rdt*(inv_thickness(i)*sty+v_rhs_ice(i)) + ice%beta_evp*v_ice_aux(i)

                    !solve (Coriolis and water stress are treated implicitly)
                    det = minval(bc_index_nod2D(edges(:,i))) / ((1.0_WP+ice%beta_evp+drag)**2 + (rdt*cor)**2)

                    u_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsu +rdt*cor*rhsv)
                    v_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsv -rdt*cor*rhsu)
                end if
            end do ! --> do i=1,myDim_edge2D
        end if
!$OMP END DO
        !_______________________________________________________________________
        ! apply sea ice velocity boundary condition
!$OMP DO
        if (ice_vplace == 0) then
            do ed=1,myDim_edge2D
                !___________________________________________________________________
                ! apply coastal sea ice velocity boundary conditions
                if (myList_edge2D(ed) > edge2D_in) then
                   do n=1, 2
#if defined(_OPENMP)
                    call omp_set_lock  (partit%plock(edges(n, ed)))
#endif
                    u_ice_aux(edges(n,ed))=0.0_WP
                    v_ice_aux(edges(n,ed))=0.0_WP
#if defined(_OPENMP)
                    call omp_unset_lock(partit%plock(edges(n,ed)))
#endif
                   end do
                end if

                !___________________________________________________________________
                ! apply sea ice velocity boundary conditions at cavity-ocean edge
                if (use_cavity) then
                    if ( (ulevels(edge_tri(1,ed))>1) .or. &
                        ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then
                        do n=1, 2
#if defined(_OPENMP)
                           call omp_set_lock  (partit%plock(edges(n, ed)))
#endif
                           u_ice_aux(edges(n,ed))=0.0_WP
                           v_ice_aux(edges(n,ed))=0.0_WP
#if defined(_OPENMP)
                           call omp_unset_lock(partit%plock(edges(n,ed)))
#endif
                        end do
                    end if
                end if
            end do ! --> do ed=1,myDim_edge2D
        else if (ice_vplace == 1) then
            do ed=1,myDim_edge2D
                !___________________________________________________________________
                ! apply coastal sea ice velocity boundary conditions
                if (myList_edge2D(ed) > edge2D_in) then
#if defined(_OPENMP)
                    call omp_set_lock  (partit%plock(ed))
#endif
                    u_ice_aux(ed)=0.0_WP
                    v_ice_aux(ed)=0.0_WP
#if defined(_OPENMP)
                    call omp_unset_lock(partit%plock(ed))
#endif
                end if

                !___________________________________________________________________
                ! apply sea ice velocity boundary conditions at cavity-ocean edge
                if (use_cavity) then
                    if (ulevels(edge_tri(1,ed))>1) then
                    
                        u_ice_aux(ed)=0.0_WP
                        v_ice_aux(ed)=0.0_WP
                    else if (edge_tri(2,ed)>0) then
                        if (ulevels(edge_tri(2,ed))>1) then

                            u_ice_aux(ed)=0.0_WP
                            v_ice_aux(ed)=0.0_WP
                        endif
                    endif

                end if
            end do ! --> do ed=1,myDim_edge2D
        end if
!$OMP END DO
        !_______________________________________________________________________
        if (ice_vplace == 0) then
!$OMP MASTER
            call exchange_nod_begin(u_ice_aux, v_ice_aux, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
            do row=1, myDim_nod2d
               u_rhs_ice(row)=0.0_WP
               v_rhs_ice(row)=0.0_WP
            end do
!$OMP END DO
!$OMP MASTER
            call exchange_nod_end(partit)
!$OMP END MASTER
!$OMP BARRIER
        else if (ice_vplace == 1) then
!$OMP MASTER
            call exchange_edge2D(u_ice_aux, partit)
            call exchange_edge2D(v_ice_aux, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
            do row=1, myDim_edge2D
                u_rhs_ice(row)=0.0_WP
                v_rhs_ice(row)=0.0_WP
            end do
!$OMP END DO
        end if
!$OMP END PARALLEL
    end do ! --> do shortstep=1, steps

!$OMP PARALLEL DO
    do row=1, myDim+eDim
       u_ice(row)=u_ice_aux(row)
       v_ice(row)=v_ice_aux(row)
    end do
!$OMP END PARALLEL DO

    ! interpolate velocities from edges to nodes if necessary
    if (ice_vplace == 0) then
        u_ice_nod = u_ice
        v_ice_nod = v_ice
    else if (ice_vplace == 1) then
        u_ice_nod = 0
        v_ice_nod = 0

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

    end if


end subroutine EVPdynamics_m
!
!
!_______________________________________________________________________________
! aEVP implementation: Similar to mEVP, but alpha is variable.
! The subroutines involved are with _a.
! EVP stability parameter alpha is computed at each element
! aEVP implementation
! SD, 13.02.2017
subroutine find_alpha_field_a(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_config
#if defined (__icepack)
    use icedrv_main,   only: strength
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: elem, elnodes(3)
    real(kind=WP)            :: dx(3), dy(3), msum, asum
    real(kind=WP)            :: eps1, eps2, pressure, delta
    real(kind=WP)            :: val3, meancos, usum, vsum, vale
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux
    real(kind=WP), dimension(:), pointer  :: alpha_evp_array
    real(kind=WP)              , pointer  :: rhoice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    eps11           => ice%work%eps11(:)
    eps12           => ice%work%eps12(:)
    eps22           => ice%work%eps22(:)
    sigma11         => ice%work%sigma11(:)
    sigma12         => ice%work%sigma12(:)
    sigma22         => ice%work%sigma22(:)
    u_ice_aux       => ice%uice_aux(:)
    v_ice_aux       => ice%vice_aux(:)
    alpha_evp_array => ice%alpha_evp_array(:)
    rhoice          => ice%thermo%rhoice

    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
    vale=1.0_WP/(ice%ellipse**2)
    do elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1) cycle

        msum=sum(m_ice(elnodes))*val3
        if(msum<=0.01_WP) cycle !DS
        asum=sum(a_ice(elnodes))*val3

        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        ! METRICS:
        vsum=sum(v_ice_aux(elnodes))
        usum=sum(u_ice_aux(elnodes))
        meancos=metric_factor(elem)
        !
        ! ====== Deformation rate tensor on element elem:
        eps11(elem)=sum(dx*u_ice_aux(elnodes))
        eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
        eps22(elem)=sum(dy*v_ice_aux(elnodes))
        eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
        eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics

        ! ======= Switch to eps1,eps2
        eps1=eps11(elem)+eps22(elem)
        eps2=eps11(elem)-eps22(elem)

        ! ====== moduli:
        delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
        delta=sqrt(delta)

#if defined (__icepack)
        pressure = sum(strength(elnodes))*val3/(delta+ice%delta_min)/msum
#else
        pressure = ice%pstar*exp(-ice%c_pressure*(1.0_WP-asum))/(delta+ice%delta_min) ! no multiplication
                                                                       ! with thickness (msum)
#endif
        !adjust c_aevp such, that alpha_evp_array and beta_evp_array become in acceptable range
        alpha_evp_array(elem)=max(50.0_WP,sqrt(ice%ice_dt*ice%c_aevp*pressure/rhoice/elem_area(elem)))
        ! /voltriangle(elem) for FESOM1.4
        ! We do not allow alpha to be too small!
    end do !--> do elem=1,myDim_elem2D
end subroutine find_alpha_field_a
!
!
!_______________________________________________________________________________
! Internal stress tensor
! New implementation following Boullion et al, Ocean Modelling 2013.
! and Kimmritz et al., Ocean Modelling 2016
! SD, 14.02.2017
subroutine stress_tensor_a(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_param
    use mod_mesh
    use g_config
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                   :: elem, elnodes(3)
    real(kind=WP)             :: dx(3), dy(3), msum, asum
    real(kind=WP)             :: eps1, eps2, pressure, delta
    real(kind=WP)             :: val3, meancos, usum, vsum, vale
    real(kind=WP)             :: det1, det2, r1, r2, r3, si1, si2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux
    real(kind=WP), dimension(:), pointer  :: alpha_evp_array
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    eps11           => ice%work%eps11(:)
    eps12           => ice%work%eps12(:)
    eps22           => ice%work%eps22(:)
    sigma11         => ice%work%sigma11(:)
    sigma12         => ice%work%sigma12(:)
    sigma22         => ice%work%sigma22(:)
    u_ice_aux       => ice%uice_aux(:)
    v_ice_aux       => ice%vice_aux(:)
    alpha_evp_array => ice%alpha_evp_array(:)

    !___________________________________________________________________________
    val3=1.0_WP/3.0_WP
    vale=1.0_WP/(ice%ellipse**2)
    do elem=1,myDim_elem2D
        !_______________________________________________________________________
        ! if element has any cavity node skip it
        if (ulevels(elem) > 1) cycle

        !_______________________________________________________________________
        det2=1.0_WP/(1.0_WP+alpha_evp_array(elem))     ! Take alpha from array
        det1=alpha_evp_array(elem)*det2

        elnodes=elem2D_nodes(:,elem)

        msum=sum(m_ice(elnodes))*val3
        if(msum<=0.01_WP) cycle !DS
        asum=sum(a_ice(elnodes))*val3

        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        ! METRICS:
        vsum=sum(v_ice_aux(elnodes))
        usum=sum(u_ice_aux(elnodes))
        meancos=metric_factor(elem)
        !
        ! ====== Deformation rate tensor on element elem:
        eps11(elem)=sum(dx*u_ice_aux(elnodes))
        eps11(elem)=eps11(elem)-val3*vsum*meancos                !metrics
        eps22(elem)=sum(dy*v_ice_aux(elnodes))
        eps12(elem)=0.5_WP*sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes))
        eps12(elem)=eps12(elem)+0.5_WP*val3*usum*meancos          !metrics

        ! ======= Switch to eps1,eps2
        eps1=eps11(elem)+eps22(elem)
        eps2=eps11(elem)-eps22(elem)

        ! ====== moduli:
        delta=eps1**2+vale*(eps2**2+4.0_WP*eps12(elem)**2)
        delta=sqrt(delta)

#if defined (__icepack)
        pressure = sum(strength(elnodes))*val3/(delta+ice%delta_min)
#else
        pressure=ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))/(delta+ice%delta_min)
#endif

        r1=pressure*(eps1-delta)
        r2=pressure*eps2*vale
        r3=pressure*eps12(elem)*vale
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*si1+det2*r1
        si2=det1*si2+det2*r2
        sigma12(elem)=det1*sigma12(elem)+det2*r3
        sigma11(elem)=0.5_WP*(si1+si2)
        sigma22(elem)=0.5_WP*(si1-si2)

#if defined (__icepack)
        rdg_conv_elem(elem)  = -min((eps11(elem)+eps22(elem)),0.0_WP)
        rdg_shear_elem(elem) = 0.5_WP*(delta - abs(eps11(elem)+eps22(elem)))
#endif
    end do ! --> do elem=1,myDim_elem2D
    ! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of
    ! Boullion et al Ocean Modelling 2013, but in an implicit mode:
    ! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
    ! and similarly for si2 and sigma12
end subroutine stress_tensor_a
!
!
!_______________________________________________________________________________
! assemble rhs and solve for ice velocity
! New implementation based on Bouillion et al. Ocean Modelling 2013
! and Kimmritz et al., Ocean Modelling  2016
! SD 14.02.17
subroutine EVPdynamics_a(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    USE o_arrays
    use o_PARAM
    use g_config, only: use_cavity
    use g_comm_auto
    use ice_maEVP_interfaces
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif
    implicit none
    type(t_ice),    intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer          :: steps, shortstep, i, ed, n
    real(kind=WP)    :: rdt, drag, det, fc
    real(kind=WP)    :: thickness, inv_thickness, umod, rhsu, rhsv
    REAL(kind=WP)    :: t0,t1, t2, t3, t4, t5, t00, txx
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: u_rhs_ice, v_rhs_ice
    real(kind=WP), dimension(:), pointer  :: u_w, v_w
    real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
    real(kind=WP), dimension(:), pointer  :: u_ice_aux, v_ice_aux
    real(kind=WP), dimension(:), pointer  :: beta_evp_array
    real(kind=WP), dimension(:), pointer  :: ice_strength
    real(kind=WP)              , pointer  :: rhoice, rhosno
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice(:)
    v_ice           => ice%vice(:)
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    m_snow          => ice%data(3)%values(:)
    u_rhs_ice       => ice%uice_rhs(:)
    v_rhs_ice       => ice%vice_rhs(:)
    u_w             => ice%srfoce_u(:)
    v_w             => ice%srfoce_v(:)
    stress_atmice_x => ice%stress_atmice_x
    stress_atmice_y => ice%stress_atmice_y
    u_ice_aux       => ice%uice_aux(:)
    v_ice_aux       => ice%vice_aux(:)
    beta_evp_array  => ice%beta_evp_array(:)
    ice_strength    => ice%work%ice_strength(:)
    rhoice          => ice%thermo%rhoice
    rhosno          => ice%thermo%rhosno

    !___________________________________________________________________________
    steps=ice%evp_rheol_steps
    rdt=ice%ice_dt
    u_ice_aux=u_ice    ! Initialize solver variables
    v_ice_aux=v_ice
    call ssh2rhs(ice, partit, mesh)

    !___________________________________________________________________________
    ! Populate the nodal ice_strength diagnostic. With icepack, use the per-
    ! node icepack strength field directly; otherwise the canonical Hibler
    ! (1979) P = P*·h·exp(-C(1-A)) from per-node m_ice and a_ice. Exposed
    ! as the 'strength_ice' output stream. aEVP evaluates the per-element
    ! strength inline in its iteration.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n)
    do n = 1, myDim_nod2D
        ice_strength(n) = 0.0_WP
        if (ulevels_nod2D(n) > 1) cycle
        if (m_ice(n) <= 0._WP .or. a_ice(n) <= 0._WP) cycle
#if defined (__icepack)
        ice_strength(n) = strength(n)
#else
        ice_strength(n) = ice%pstar*m_ice(n)*exp(-ice%c_pressure*(1.0_WP-a_ice(n)))
#endif
    end do
!$OMP END PARALLEL DO

#if defined (__icepack)
    rdg_conv_elem(:)  = 0.0_WP
    rdg_shear_elem(:) = 0.0_WP
#endif

    do shortstep=1, steps
        call stress_tensor_a(ice, partit, mesh)
        call stress2rhs_m(ice, partit, mesh)    ! _m=_a, so no _m version is the only one!
        do i=1,myDim_nod2D

            !___________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels_nod2d(i)>1) cycle

            thickness=(rhoice*m_ice(i)+rhosno*m_snow(i))/max(a_ice(i),0.01_WP)
            thickness=max(thickness, 9.0_WP)   ! Limit if it is too small (0.01 m)
            inv_thickness=1.0_WP/thickness

            umod=sqrt((u_ice_aux(i)-u_w(i))**2+(v_ice_aux(i)-v_w(i))**2)
            drag=rdt*ice%cd_oce_ice*umod*density_0*inv_thickness

            !rhs for water stress, air stress, and u_rhs_ice/v (internal stress + ssh)
            rhsu=u_ice(i)+drag*u_w(i)+rdt*(inv_thickness*stress_atmice_x(i)+u_rhs_ice(i))
            rhsv=v_ice(i)+drag*v_w(i)+rdt*(inv_thickness*stress_atmice_y(i)+v_rhs_ice(i))

            rhsu=beta_evp_array(i)*u_ice_aux(i)+rhsu
            rhsv=beta_evp_array(i)*v_ice_aux(i)+rhsv
            !solve (Coriolis and water stress are treated implicitly)
            fc=rdt*mesh%coriolis_node(i)
            det=(1.0_WP+beta_evp_array(i)+drag)**2+fc**2
            det=bc_index_nod2D(i)/det
            u_ice_aux(i)=det*((1.0_WP+beta_evp_array(i)+drag)*rhsu+fc*rhsv)
            v_ice_aux(i)=det*((1.0_WP+beta_evp_array(i)+drag)*rhsv-fc*rhsu)
        end do

        !_______________________________________________________________________
        ! apply sea ice velocity boundary condition
        do ed=1,myDim_edge2D
            !___________________________________________________________________
            ! apply coastal sea ice velocity boundary conditions
            if(myList_edge2D(ed) > edge2D_in) then
                u_ice_aux(edges(:,ed))=0.0_WP
                v_ice_aux(edges(:,ed))=0.0_WP
            end if

            !___________________________________________________________________
            ! apply sea ice velocity boundary conditions at cavity-ocean edge
            if (use_cavity) then
                if ( (ulevels(edge_tri(1,ed))>1) .or. &
                    ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then
                    u_ice_aux(edges(1:2,ed))=0.0_WP
                    v_ice_aux(edges(1:2,ed))=0.0_WP
                end if
            end if
        end do ! --> do ed=1,myDim_edge2D

        call exchange_nod(u_ice_aux, v_ice_aux, partit)
    end do

    u_ice=u_ice_aux
    v_ice=v_ice_aux

    call find_alpha_field_a(ice, partit, mesh)             ! alpha_evp_array is initialized with alpha_evp;
                                        ! At this stage we already have non-trivial velocities.
    call find_beta_field_a(ice, partit, mesh)
end subroutine EVPdynamics_a
!
!
!_______________________________________________________________________________
! beta_evp_array is defined at nodes, and this is the only
! reason we need it in addition to alpha_evp_array (we work with
! alpha=beta, and keep different names for generality; mEVP can work with
! alpha \ne beta, but not aEVP).
subroutine find_beta_field_a(ice, partit, mesh)
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE MOD_ICE
    use o_param
    Implicit none
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_ice)   , intent(inout), target :: ice
    !___________________________________________________________________________
    integer :: n
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: alpha_evp_array, beta_evp_array
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    alpha_evp_array => ice%alpha_evp_array(:)
    beta_evp_array  => ice%beta_evp_array(:)

    !___________________________________________________________________________
    DO n=1, myDim_nod2D
       !________________________________________________________________________
       ! if element has any cavity node skip it
       if (ulevels_nod2d(n)>1) cycle

       ! ==============
       ! FESOM1.4 and stand-alone FESIM
       ! beta_evp_array(n) =  maxval(alpha_evp_array(nod_in_elem2D(n)%addresses(1:nod_in_elem2D(n)%nmb)))
       ! ==============
       ! FESOM2.0
       beta_evp_array(n) =  maxval(alpha_evp_array(nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
    END DO
end subroutine find_beta_field_a
!
! ================================================================
!
