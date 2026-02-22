"""
HQIV Cosmology — Horizon-Quantized Informational Vacuum Perturbation Solver.

A pure-Python (NumPy + SciPy) Fourier-space integrator for the HQIV framework
that evolves scalar and vector perturbation equations from deep radiation era
through recombination to today.

Paper Reference: "Horizon-Quantized Informational Vacuum (HQIV): 
A Covariant Baryon-Only Cosmological Framework from Quantised Inertia"

Main Components
---------------
background : Background cosmology module
    CosmologicalBackground class for H(a), η(a), G_eff(a), β(a)
    
scalars : Scalar perturbation module
    ScalarPerturbations class for δ, θ, Φ evolution
    Implements paper eqs. 9-11 with modified inertia
    
vectors : Vector perturbation module
    VorticityPerturbations class for ω evolution
    Implements paper eq. 12 with horizon coupling
    
observables : Observables module
    CMBPowerSpectrum for C_ℓ^TT computation
    Matter power spectrum and σ₈
    
hqiv_solver : Main solver
    HQIVPerturbations master class

Quick Start
-----------
>>> from hqiv_solver import HQIVPerturbations
>>> solver = HQIVPerturbations(H0=73.2, Om_m=0.031, hqiv_on=True)
>>> solver.run()
>>> solver.plot_results()

Author: HQIV Team
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "HQIV Team"

# Import main classes
from .background import CosmologicalBackground
from .scalars import ScalarPerturbations, compute_transfer_function
from .vectors import VorticityPerturbations, VorticityScalarCoupling
from .observables import CMBPowerSpectrum, simplified_CMB_spectrum
from .hqiv_solver import HQIVPerturbations, quick_test

# Import constants
from .background import c, Mpc_m, G0_SI, Gyr_s

__all__ = [
    # Main classes
    'HQIVPerturbations',
    'CosmologicalBackground',
    'ScalarPerturbations',
    'VorticityPerturbations',
    'VorticityScalarCoupling',
    'CMBPowerSpectrum',
    
    # Functions
    'compute_transfer_function',
    'simplified_CMB_spectrum',
    'quick_test',
    
    # Constants
    'c',
    'Mpc_m',
    'G0_SI',
    'Gyr_s',
]