! HQIV covariant extension — Steven Ettinger + Mr 4.20
! Paste after "end function GetGratio" in results.f90.
!
! ! Theta(a) = 2c/H(a) in Mpc. H in km/s/Mpc => Theta = 2*c_km_s/(H0*H_over_H0) Mpc
! function GetTheta(this, a) result(Theta_Mpc)
!   class(CAMBdata) :: this
!   real(dl), intent(in) :: a
!   real(dl) :: Theta_Mpc
!   real(dl), parameter :: c_km_s = 299792.458_dl
!   Theta_Mpc = 0._dl
!   if (.not. this%HQIV_mode .or. this%hqiv_na < 1) return
!   Theta_Mpc = (2._dl * c_km_s) / (this%CP%H0 * GetHoverH0(this, a))
! end function GetTheta
!
! ! Inertia reduction m_i/m_g = 1 - a_min/|a_local|. a_min = beta*c*H. a_local_scale > 0 in same units.
! function GetInertiaFactor(this, a, a_local_scale) result(fac)
!   use constants, only : c
!   class(CAMBdata) :: this
!   real(dl), intent(in) :: a, a_local_scale
!   real(dl) :: fac, a_min
!   real(dl), parameter :: c_km_s = 299792.458_dl
!   fac = 1._dl
!   if (.not. this%CP%HQIV_covariant) return
!   if (a_local_scale <= 0._dl) return
!   a_min = this%CP%hqiv_beta * c_km_s * this%CP%H0 * GetHoverH0(this, a) * 1e-3_dl
!   fac = 1._dl - min(a_min / a_local_scale, 0.99_dl)
!   fac = max(fac, 0.01_dl)
! end function GetInertiaFactor
