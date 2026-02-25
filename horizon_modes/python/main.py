import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
N_steps = 6000
T_Pl   = 1.22e19      # GeV
T_QCD  = 1.8          # GeV  (lock-in window)
gamma_0 = 0.40        # fiducial γ from paper (3H² − γH = 8π G_eff ρ); sets scale for γ(T)
handedness_factor = 1.0 / (4 * np.pi)  # 1/(4π) conservative

print("=== Discrete Shell-by-Shell Mode Counting with Derived Vorticity Asymmetry ===\n")

# Log-spaced temperature grid (T decreases)
logT = np.linspace(np.log10(T_Pl), np.log10(T_QCD), N_steps)
T    = 10**logT

# Shell index k = Θ_local / l_Pl  (Θ_local ∝ 1/T² in RD)
k = (T_Pl / T)**2                     # starts at k=1

# New modes per shell — holographic area law (Brodie/Jacobson)
dN_new = 8 * np.pi * k

# === γ(T) from the vorticity term in the paper (main.tex) ===
# Paper: 3H² − γH = 8π G_eff(φ) ρ_tot. The −γH term is the horizon/inertia (vorticity) term.
# φ = 2c²/Θ_local, in RD Θ_local ∝ 1/T² ⇒ φ ∝ T². So φ/φ_Pl ∝ (T/T_Pl)².
# Define γ per step from the vorticity term: γ(φ) = γ_0 × (φ/φ_Pl) ⇒ γ(T) = γ_0 × (T/T_Pl)².
phi_ratio = (T / T_Pl)**2             # φ/φ_Pl in RD
gamma_T = gamma_0 * phi_ratio         # γ(T) from paper's vorticity term

# === DERIVATION OF ASYMMETRY AMPLITUDE FROM THE ACTION ===
# 1. ϕ = 2c² / Θ_local  →  ϕ ∝ T²
# 2. From action: ∂f/∂ϕ ≈ -c/(6 (a_loc + cϕ/6))   (Brodie pure form)
#    In the transition regime near QCD we approximate ∂f/∂ϕ ≈ -1/(6 ϕ) (normalised)
df_dphi = -1.0 / (6.0 * (1.0 + (T/T_Pl)**2))

# 3. Typical |∇ϕ| from horizon discreteness; coefficient γ(T) from vorticity term (paper).
#    |∇ϕ|/ϕ_Pl ~ γ(T) × (T/T_Pl)^4 (modified Einstein eq, γ now step-dependent).
grad_phi = gamma_T * (T / T_Pl)**4

# 4. Vorticity source strength per mode = |∂f/∂ϕ * (k × ∇ϕ)|
#    Handedness average gives factor ~1/(4π) (tunable via handedness_factor)
delta = np.abs(df_dphi * grad_phi) * handedness_factor

# Accumulate totals
cum_total_modes = np.cumsum(dN_new)
cum_net_baryons = np.cumsum(dN_new * delta)

# Final result at QCD lock-in
eta = cum_net_baryons[-1] / cum_total_modes[-1]

# Convergence check: compare with 2× steps (steps are sufficient if relative change is tiny)
def compute_eta(n_steps):
    logT_ = np.linspace(np.log10(T_Pl), np.log10(T_QCD), n_steps)
    T_ = 10**logT_
    k_ = (T_Pl / T_)**2
    dN_ = 8 * np.pi * k_
    gamma_T_ = gamma_0 * (T_ / T_Pl)**2   # γ(T) from vorticity term
    df_ = -1.0 / (6.0 * (1.0 + (T_/T_Pl)**2))
    g_ = gamma_T_ * (T_ / T_Pl)**4
    d_ = np.abs(df_ * g_) * handedness_factor
    return np.sum(dN_ * d_) / np.sum(dN_)

eta_double = compute_eta(2 * N_steps)
rel_change = abs(eta_double - eta) / (eta + 1e-100)
steps_ok = rel_change < 0.01

print(f"Steps                  : {N_steps}")
print(f"Convergence (2× steps) : {rel_change:.2%} change in η → {'OK' if steps_ok else 'try more steps'}")
print(f"Final total modes      : {cum_total_modes[-1]:.3e}")
print(f"Final net baryons      : {cum_net_baryons[-1]:.3e}")
print(f"Derived η at QCD       : {eta:.3e}")
print(f"Observed target        : 6.1e-10")
print(f"Factor vs observed     : {eta / 6.1e-10:.2e}x")
if eta < 6.1e-10:
    print("\n--- Suggestions to approach target ---")
    print("  • γ(T) = γ_0 (T/T_Pl)² from paper vorticity term; try larger γ_0 or different γ(φ) form.")
    print("  • Try handedness_factor = 1/(2π) or 1/π if vorticity is better aligned.")
    print("  • Check ∂f/∂ϕ normalisation and coupling strength near QCD.")
    print("  • Confirm η definition: baryon-to-photon may need mode-to-photon conversion.")

# Plot: η(T) and γ(T) from paper vorticity term
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
ax1.semilogx(T, cum_net_baryons / cum_total_modes, 'b-', lw=2.5, label='Cumulative η from vorticity')
ax1.axvline(T_QCD, color='red', ls='--', label='QCD lock-in window')
ax1.set_ylabel('Net baryon-to-photon ratio η')
ax1.set_title('Emergent Baryon Asymmetry from Pre-Plasma Horizon Vorticity')
ax1.grid(True, which='both')
ax1.legend()
ax2.semilogx(T, gamma_T, 'g-', lw=1.5, label=r'γ(T) = γ₀ (T/T_Pl)² from vorticity term (paper)')
ax2.axvline(T_QCD, color='red', ls='--', alpha=0.7)
ax2.set_xlabel('Temperature [GeV]')
ax2.set_ylabel('γ(T)')
ax2.grid(True, which='both')
ax2.legend()
plt.tight_layout()
if plt.get_backend().lower() != 'agg':
    plt.show()
else:
    plt.savefig('eta_vs_T.png', dpi=120)
    plt.close()