# HQIV: add to Python camb/model.py so the CAMBparams ctypes struct matches Fortran.
# Add these entries to CAMBparams._fields_ in the SAME order as in fortran/model.f90
# (immediately after "Alens"):
#
#         ("Alens", c_double, "non-physical scaling amplitude for the CMB lensing spectrum power"),
#         ("HQIV", c_bool, "Use HQIV cosmology (tabulated H(a), baryons only, horizon cutoff)"),
#         ("hqiv_Ha_file", c_char * 1024, "Path to table: a, H_over_H0 (Ini_max_string_len=1024)"),
#         ("hqiv_cs2_fac", c_double, "Effective sound speed factor placeholder (~0.95-1.05)"),
#         ("hqiv_l_cut", c_double, "Low-ell cutoff scale for exp(-(l/l_cut)^1.8)"),
#         ("hqiv_beta", c_double, "Horizon motive beta (~0.78)"),
#         ("MassiveNuMethod", c_int, ...),
#
# If your camb package uses a different constant for string length, use that instead of 1024.
