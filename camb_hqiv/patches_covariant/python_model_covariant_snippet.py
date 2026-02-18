# HQIV covariant extension — Steven Ettinger + Mr 4.20
# Add to camb/model.py in CAMBparams._fields_, in the SAME order as fortran/model.f90 (after hqiv_beta):
#
#         ("hqiv_beta", c_double, "Horizon motive beta (~0.78)"),
#         ("HQIV_covariant", c_bool, "Full covariant HQIV: horizon metric, inertia reduction, varying G, horizon cutoff"),
#         ("MassiveNuMethod", c_int, ...),
#
