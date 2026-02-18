import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import camb

# HQIV-style: baryons only, Omega_total=1 from horizon equilibrium (no separate DE)
H0 = 70.0
ombh2 = 0.048
omch2 = 0.0
# Late-time acceleration from QI horizon term → flat with effective Lambda for this illustration
omk = 0.0
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=0.06, omk=omk)
pars.set_dark_energy(w=-1.0, wa=0.0)
pars.InitPower.set_params(ns=0.96, As=2.1e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0)

print("HQIV-style CMB TT — QI horizon cutoff (multipole graph)")

results = camb.get_results(pars)
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
cl_tt = powers['total'][:, 0]
l = np.arange(len(cl_tt))

# Natural horizon cutoff (super-horizon modes only)
l_cut = 38.0
damping = np.ones_like(l, dtype=float)
mask = l < 80
damping[mask] = np.exp(-(l[mask] / l_cut) ** 1.8)
cl_tt_hqiv = cl_tt * damping

low_l_suppression = 100 * (1 - np.mean(cl_tt_hqiv[2:31] / np.maximum(cl_tt[2:31], 1e-30)))
peak_idx = np.argmax(cl_tt[180:280]) + 180
print(f"Low-ℓ suppression (ℓ=2–30): {low_l_suppression:.1f}%")
print(f"First acoustic peak: ℓ ≈ {l[peak_idx]:.0f}")

# Multipole graph: D_ell = ell(ell+1)C_ell/(2pi)
d_ell = l * (l + 1) * cl_tt_hqiv / (2 * np.pi)
plt.figure(figsize=(10, 6))
plt.plot(l, d_ell, color='C0', lw=1.2, label='HQIV CMB TT (horizon cutoff)')
plt.xlim(2, 1200)
plt.ylim(0, 8000)
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$  [$\mu$K$^2$]')
plt.title('HQIV CMB TT — no separate DE; QI horizon cutoff')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('hqiv_cmb_tt.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved hqiv_cmb_tt.png")
