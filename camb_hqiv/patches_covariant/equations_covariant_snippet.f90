! Full covariant HQIV — equations.f90 (GaugeInterface / derivs).
! NOTE: Horizon cutoff exp(-(k/k_cut)^1.8) is already implemented in cmbmain.f90 CalcScalCls.
!
! 1) Effective sound speed (baryon-photon fluid): c_s,eff^2 = c_s^2 / (1 - beta*c*H/|a_fluid|)
!    Where cs2 is used (e.g. in derivs, or in Thermo_values call), after getting cs2:
!    if (CP%HQIV_covariant .and. State%HQIV_mode) then
!      Theta_a = GetTheta(State, a)
!      a_fluid = max(cs2*k/a, State%CP%H0*GetHoverH0(State,a)*1e-3_dl*299792.458_dl)
!      a_min = State%CP%hqiv_beta * 299792.458_dl * State%CP%H0 * GetHoverH0(State,a) * 1e-3_dl
!      denom = 1._dl - a_min / max(a_fluid, 1e-30_dl)
!      cs2 = cs2 / max(denom, 0.01_dl)
!    end if
!
! 2) Inertia factor in continuity / Euler
!    Where the baryon velocity or density evolution uses inertial mass (e.g. vbdot, delta_b),
!    multiply the relevant inertial terms by GetInertiaFactor(State, a, a_local_scale).
!    a_local_scale = max(cs*k/a, H_phys) or similar in consistent units.
!
! 3) Poisson equation (varying G)
!    Where the metric potential is sourced by 8*pi*G*delta_rho, multiply the source by
!    GetGratio(State, a) when CP%HQIV or CP%HQIV_covariant.
!
! 4) Horizon cutoff in CMB source (SW + ISW)
!    In the source function that feeds the C_l integral (e.g. in the visibility or
!    time-integrated source), for each k and tau: k_cut = 2*pi/GetTheta(State,a).
!    Apply damping: if (k < k_cut) source = source * exp(-(k/k_cut)**1.8)
!    Use a = 1/(1+z) from tau via State%tau0 and the time step.
