! ========== HQIV: paste into results.f90 ==========
! 1) In type CAMBdata, after "real(dl) ThermoDerivedParams(nthermo_derived)" add:
!
!    ! HQIV: tabulated H(a)/H0
!    logical :: HQIV_mode = .false.
!    integer :: hqiv_na = 0
!    real(dl), allocatable :: hqiv_a(:), hqiv_H_over_H0(:)
!
!    In the same type's "procedure" list, add:  procedure :: LoadHQIVTable
!
! 2) Add these subroutines in the "contains" section of module results (e.g. before CAMBdata_SetParams):
!
! subroutine LoadHQIVTable(this)
!   use constants, only : c
!   class(CAMBdata) :: this
!   integer :: i, ierr, na
!   character(LEN=512) :: line
!   real(dl) :: a, Hov
!   if (.not. this%CP%HQIV) return
!   if (allocated(this%hqiv_a)) deallocate(this%hqiv_a, this%hqiv_H_over_H0)
!   open(unit=99, file=trim(adjustl(this%CP%hqiv_Ha_file)), status='old', action='read', iostat=ierr)
!   if (ierr /= 0) then
!     call GlobalError('HQIV: cannot open ' // trim(adjustl(this%CP%hqiv_Ha_file)), error_unsupported_params)
!     return
!   end if
!   read(99, *, iostat=ierr) line   ! skip header
!   na = 0
!   do
!     read(99, *, iostat=ierr) a, Hov
!     if (ierr /= 0) exit
!     na = na + 1
!   end do
!   rewind(99)
!   read(99, *, iostat=ierr) line
!   allocate(this%hqiv_a(na), this%hqiv_H_over_H0(na))
!   do i = 1, na
!     read(99, *) this%hqiv_a(i), this%hqiv_H_over_H0(i)
!   end do
!   close(99)
!   this%hqiv_na = na
!   this%HQIV_mode = .true.
! end subroutine LoadHQIVTable
!
! function GetHoverH0(this, a) result(H_over_H0)
!   class(CAMBdata) :: this
!   real(dl), intent(in) :: a
!   real(dl) :: H_over_H0
!   integer :: j
!   real(dl) :: loga, loga1, loga2, w
!   H_over_H0 = 1._dl
!   if (.not. this%HQIV_mode .or. this%hqiv_na < 2) return
!   if (a <= this%hqiv_a(1)) then
!     H_over_H0 = this%hqiv_H_over_H0(1)
!     return
!   end if
!   if (a >= this%hqiv_a(this%hqiv_na)) then
!     H_over_H0 = this%hqiv_H_over_H0(this%hqiv_na)
!     return
!   end if
!   do j = 1, this%hqiv_na - 1
!     if (a >= this%hqiv_a(j) .and. a <= this%hqiv_a(j+1)) then
!       loga1 = log(max(this%hqiv_a(j), 1.d-30))
!       loga2 = log(this%hqiv_a(j+1))
!       loga = log(a)
!       w = (loga - loga1) / (loga2 - loga1)
!       H_over_H0 = (1._dl - w) * this%hqiv_H_over_H0(j) + w * this%hqiv_H_over_H0(j+1)
!       return
!     end if
!   end do
! end function GetHoverH0
!
! function GetGratio(this, a) result(G_ratio)
!   ! G(a)/G0 = (Theta0/Theta(a))^0.6 = (H(a)/H0)^0.6
!   class(CAMBdata) :: this
!   real(dl), intent(in) :: a
!   real(dl) :: G_ratio
!   G_ratio = GetHoverH0(this, a) ** 0.6_dl
! end function GetGratio
!
! 3) In CAMBdata_SetParams, after "this%grhok=this%grhocrit*this%CP%omk" add:
!
!   if (this%CP%HQIV) then
!     this%CP%omch2 = 0._dl
!     this%grhoc = 0._dl
!     call this%LoadHQIVTable()
!     if (this%hqiv_na > 0) then
!       this%Omega_de = 1._dl - (this%CP%ombh2 + this%CP%omnuh2)/h2 - this%CP%omk &
!         - (this%grhornomass + this%grhog)/this%grhocrit
!       this%grhov = this%grhocrit * this%Omega_de
!     end if
!   end if
!
! 4) In CAMBdata_Hofz, at the start of the function body (after "a = 1/(1+z)"), add:
!
!   if (this%HQIV_mode .and. this%hqiv_na > 0) then
!     CAMBdata_Hofz = (this%CP%H0/299792.458_dl) * GetHoverH0(this, a)
!     return
!   end if
!
! 5) Low-ell damping: apply exp(-(ell/ell_cut)^1.8) to Cl for ell < 80 when HQIV.
!    Find where Cl_scalar(l, C_Temp) is set (e.g. in subroutines that fill CLdata%Cl_scalar).
!    After the normal assignment, add (pseudocode):
!      if (CP%HQIV .and. l < 80 .and. l >= 2) then
!        fac = exp(-(real(l,dl)/CP%hqiv_l_cut)**1.8)
!        Cl_scalar(l, :) = Cl_scalar(l, :) * fac
!      end if
