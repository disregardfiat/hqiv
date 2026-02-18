! HQIV covariant extension — Steven Ettinger + Mr 4.20
! Horizon cutoff for Sachs-Wolfe + ISW: apply to C_l source integrand when HQIV_covariant = T.
!
! In the routine that computes the CMB source (e.g. in equations.f90 or the module that
! fills the time/k source array for the line-of-sight integral), for each (k, tau):
!
!   if (CP%HQIV_covariant .and. State%hqiv_na > 0) then
!     Theta_Mpc = GetTheta(State, a)
!     if (Theta_Mpc > 0._dl) then
!       k_cut = 2._dl * const_pi / Theta_Mpc
!       if (k < k_cut) source_k_tau = source_k_tau * exp(-(k/max(k_cut,1e-30_dl))**1.8)
!     end if
!   end if
!
! a is obtained from tau via the background (e.g. from the time step table).
