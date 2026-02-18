! ========== HQIV: paste into camb.f90 in CAMB_ReadParams ==========
!
! After reading "P%h0 = Ini%Read_Double('hubble')" (or after the use_physical block), add:
!
!   P%HQIV = Ini%Read_Logical('HQIV', .false.)
!   if (P%HQIV) then
!     call Ini%Read('hqiv_Ha_file', P%hqiv_Ha_file)
!     call Ini%Read('hqiv_cs2_fac', P%hqiv_cs2_fac)
!     call Ini%Read('hqiv_l_cut', P%hqiv_l_cut)
!     call Ini%Read('hqiv_beta', P%hqiv_beta)
!   end if
!
! If your CAMB version uses different Ini%Read APIs (e.g. Read_String_Default for the file),
! adjust accordingly. For Read_String_Default: P%hqiv_Ha_file = Ini%Read_String_Default('hqiv_Ha_file', 'hqiv_Ha.txt').
