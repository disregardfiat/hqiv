! ========== HQIV: paste into equations.f90 ==========
!
! At the top of function dtauda (after "real(dl) :: dtauda, grhoa2, grhov_t"), add:
!
!   use constants, only : c
!
! Then at the very start of the executable part of dtauda (before "call this%CP%DarkEnergy%BackgroundDensityAndPressure"),
! insert:
!
!   if (this%CP%HQIV .and. this%hqiv_na > 0) then
!     ! dtau/da = 1/(a^2 H) with H in 1/Mpc. H0 in km/s/Mpc, so H0 in 1/Mpc = H0/(c in km/s)
!     dtauda = (c/1000._dl) / (a**2 * this%CP%H0 * GetHoverH0(this, a))
!     return
!   end if
!
! Ensure "use results" is present so GetHoverH0 (from results) is visible. If dtauda is outside the
! module that uses results, you may need to add an interface or make GetHoverH0 a module procedure
! in results and use results in the scope where dtauda is defined.
