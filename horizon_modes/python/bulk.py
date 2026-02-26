# forward_4d_hqiv.py  — Baryogenesis to lock-in, then modified Friedmann to T_cmb
#
# Phase 1: Run discrete baryogenesis only up to the QCD lock-in (T ~ 1.8 GeV).
#          Output: η, Ω_k^true, and boundary state (a_lock, rho_r, rho_b, T_lock).
# Phase 2: Integrate the correct HQIV modified Friedmann (3H² - γ_eff H = 8π G_eff ρ - 3K/a²)
#          with curvature K = -Ω_k^true H0² from lock-in down to T_cmb. Write table for CLASS.
#
# CLASS then reads this table, sets a=1 at the row where T = T_cmb, and runs thermo + perturbations.
import numpy as np
from discrete_baryogenesis_horizon import (
    run_simulation,
    T_PL_GEV,
    T_QCD_GEV,
    GAMMA,
    curvature_imprint_energy,
)

# CMB temperature today (K) — CLASS runs to this; table must include this row
T_CMB_K = 2.725
# Curvature: paper value from discrete null-lattice. Same mechanism sources both
# Ω_k^true and the per-shell δE(m) → η; in the code, run_simulation() uses this
# value to scale the curvature imprint (mild circularity: η is generated with
# target Ω_k as input). For first-principles η only, call run_simulation(…,
# use_omega_k_amplitude=False); then Ω_k would need to co-emerge from a shell
# integral (see discrete_baryogenesis_horizon.py).
OMEGA_K_TRUE = 0.0098
# Recombination scale factor for γ_eff(a) ramp
A_REC = 1.0 / 1090.0
DA_REC = 0.02 * A_REC
# GeV → Kelvin for radiation temperature
GEV_TO_K = 1.160e13
# Convert comoving density in GeV^4 to CLASS units (1/Mpc)^2) so that H = sqrt(3*rho/(3-gamma)) gives 1/Mpc
# 8πG ρ = 3 H² (CLASS convention). ρ_CLASS = ρ_GeV4 * (8πG/3) in (1/Mpc)²: 1 Mpc = 1.564e38 GeV⁻¹ → 1/Mpc = 6.4e-39 GeV.
# H_GeV = sqrt(8π ρ_GeV4 / (3 M_Pl^2)); ρ_CLASS = H² so ρ_CLASS = 8π ρ_GeV4 / (3 M_Pl^2) in GeV², then (1/Mpc)² = 4.1e-77 GeV² → ρ_CLASS = ρ_GeV4 * 8π/(3*M_Pl^2) * 2.44e76 ≈ ρ_GeV4 * 2.23e15
GEV4_TO_MPC2 = 2.23e15
M_PL_GEV = 1.22e19


def _gamma_eff(a: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """γ_eff(a): 0 before recombination, γ after (smooth tanh)."""
    return gamma * 0.5 * (1.0 + np.tanh((a - A_REC) / DA_REC))


def _lockin_index(data: dict) -> int:
    """Index where lock-in is sharpest (max |d p_lock / d ln T|)."""
    p = data["p_lock"]
    lnT = np.log(data["T"] + 1e-300)
    dpl = np.gradient(p, lnT)
    return int(np.argmax(np.abs(dpl)))


def _friedmann_hqiv(
    a: float,
    rho_r_comov: float,
    rho_b_comov: float,
    gamma: float,
    H0: float,
    K: float,
) -> float:
    """
    Solve 3 H² - γ_eff H = 8πG ρ_phys - 3K/a² in natural units (G = 1/M_Pl²).
    ρ in GeV⁴, H and K in GeV² (1/time²). Returns H in GeV.
    """
    rho_r = rho_r_comov / (a ** 4)
    rho_b = rho_b_comov / (a ** 3)
    rho_phys_gev4 = rho_r + rho_b
    # 3 H² - γ_eff H0 H = 8π ρ / M_Pl² - 3K/a²  => 3 H² - γ_eff H0 H - (8π ρ/M_Pl² - 3K/a²) = 0
    source = 8.0 * np.pi * rho_phys_gev4 / (M_PL_GEV ** 2) - 3.0 * K / (a * a)
    gamma_eff = _gamma_eff(np.array([a]), gamma)[0]
    gamma_class = gamma_eff * H0
    disc = gamma_class ** 2 + 12.0 * source
    if disc < 0:
        disc = 0.0
    return (gamma_class + np.sqrt(disc)) / 6.0


def forward_4d_evolution(
    n_steps: int = 6000,
    gamma: float = GAMMA,
    T_cmb_k: float = T_CMB_K,
    omega_k_true: float = OMEGA_K_TRUE,
    seed: int = 42,
    n_loga: int = 4000,
    outpath: str = "hqiv_lattice_table.dat",
):
    """
    Run baryogenesis to lock-in, then integrate modified Friedmann to T_cmb.
    Writes (loga, rho_r_comov, rho_b_comov, T_K) in CLASS units for hqiv_emergent.
    """
    # --- Phase 1: Baryogenesis up to lock-in ---
    data = run_simulation(n_steps=n_steps, seed=seed, E_0_factor=1.0, transition_m=500)
    eta = data["eta"]
    # Use |η| for density so that ρ_b > 0 (paper predicts η > 0; sign is chiral)
    eta_for_rho = np.abs(eta)
    T_gev = data["T"]
    m = data["m"]
    R_h = data["R_h"]
    dN_new = data["dN_new"]
    bias_total = data["bias_total"]

    idx_lock = _lockin_index(data)
    T_lock_gev = float(T_gev[idx_lock])
    T_lock_k = T_lock_gev * GEV_TO_K

    cum_modes = np.cumsum(dN_new)
    cum_baryons = np.cumsum(bias_total * dN_new)
    # Comoving radiation from entropy: S ∝ cum_modes, T ∝ S^(1/3) in RD, so rho_r ∝ T^4, rho_r_comov = rho_r * a^4 ∝ T^4 * (T_ref/T)^3 = T * T_ref^3 (no) — in RD a ∝ 1/T, so rho_r_comov = rho_r * a^4 = (π²/15) T^4 * (T_ref/T)^4 = (π²/15) T_ref^4 constant. So at lock-in, rho_r_comov = (π²/15) T_lock^4 (in GeV^4) and rho_b_comov from baryon count.
    rho_r_lock_gev4 = (np.pi ** 2 / 15.0) * (T_lock_gev ** 4)
    # Baryon mass density at lock-in: n_b = η n_γ, ρ_b = η n_γ m_b. n_γ = (2ζ(3)/π²)T³. Comoving ρ_b = ρ_b a³; in RD a ∝ 1/T so ρ_b_comov ∝ η T³ a³ = η T_ref³ at lock-in.
    zeta3 = float(np.sum(1.0 / np.arange(1, 200) ** 3))
    n_gamma_lock = (2.0 * zeta3 / np.pi ** 2) * (T_lock_gev ** 3)
    m_b_gev = 0.938
    rho_b_lock_gev4 = eta_for_rho * n_gamma_lock * m_b_gev
    # a at lock-in: T ∝ 1/a => a_lock = T_cmb / T_lock (both in K).
    a_lock = T_cmb_k / T_lock_k

    # --- Phase 2: Modified Friedmann from lock-in to T_cmb ---
    # Comoving densities are constant from lock-in to today.
    rho_r_comov_gev4 = rho_r_lock_gev4 * (a_lock ** 4)
    rho_b_comov_gev4 = rho_b_lock_gev4 * (a_lock ** 3)

    # H0 closure: iterate so that at a=1, H(1) = H0. K = -Ω_k_true * H0². H in GeV from 3H² - γH = 8πρ/M_Pl².
    rho_today_gev4 = rho_r_comov_gev4 + rho_b_comov_gev4
    source_today = 8.0 * np.pi * max(1e-60, rho_today_gev4) / (M_PL_GEV ** 2)
    H0_guess = (gamma + np.sqrt(gamma ** 2 + 12.0 * source_today)) / 6.0
    for _ in range(20):
        K_guess = -omega_k_true * (H0_guess ** 2)
        H_at_1 = _friedmann_hqiv(
            1.0, rho_r_comov_gev4, rho_b_comov_gev4, gamma, H0_guess, K_guess
        )
        if abs(H_at_1 - H0_guess) < 1e-12 * H0_guess:
            break
        H0_guess = H_at_1
    H0_gev = H0_guess
    K_gev = -omega_k_true * (H0_gev ** 2)

    # Build log(a) grid from a_lock to 1
    loga_min = np.log(a_lock)
    loga_max = 0.0
    loga_arr = np.linspace(loga_min, loga_max, n_loga)
    a_arr = np.exp(loga_arr)
    T_arr_k = T_cmb_k / a_arr

    H_arr = np.zeros_like(a_arr)
    for i, a in enumerate(a_arr):
        H_arr[i] = _friedmann_hqiv(
            a, rho_r_comov_gev4, rho_b_comov_gev4, gamma, H0_gev, K_gev
        )

    # Convert to CLASS units: densities in (1/Mpc)^2 so that H = sqrt(3*rho/(3-γ)) gives 1/Mpc
    # Table stores comoving: rho_r_table = rho_r_comov (physical rho_r * a^4), rho_b_table = rho_b_comov.
    # At row i, physical rho_r = rho_r_table[i] / a^4, so table = rho_comov in CLASS units.
    # rho_comov in GeV^4 * GEV4_TO_MPC2 = rho_comov in (1/Mpc)^2.
    rho_r_comov_class = rho_r_comov_gev4 * GEV4_TO_MPC2
    rho_b_comov_class = rho_b_comov_gev4 * GEV4_TO_MPC2

    np.savetxt(
        outpath,
        np.column_stack((
            loga_arr,
            np.full_like(loga_arr, rho_r_comov_class),
            np.full_like(loga_arr, rho_b_comov_class),
            T_arr_k,
        )),
        header="loga rho_r_comov rho_b_comov T_K",
        comments="#",
    )
    print(f"Lattice table saved to {outpath} — ready for CLASS (T_cmb = {T_cmb_k} K)")

    # What CLASS will compute from the table: H0 = sqrt(3*rho_today/(3-gamma)) in 1/Mpc
    rho_today_class = rho_r_comov_class + rho_b_comov_class
    H0_class_mpc = np.sqrt(3.0 * rho_today_class / (3.0 - gamma))
    H0_km_s_Mpc = H0_class_mpc * 2.998e5  # c in km/s
    omega_m_class = rho_b_comov_class / (rho_today_class + 1e-300)
    omega_r_class = rho_r_comov_class / (rho_today_class + 1e-300)

    print("\n=== HQIV: Baryogenesis to lock-in, then Friedmann to T_cmb ===")
    print(f"Lock-in: T = {T_lock_gev:.4g} GeV ({T_lock_k:.2e} K), a_lock = {a_lock:.4e}")
    print(f"η (from baryogenesis) = {eta:.4e}")
    print(f"Ω_k^true = {omega_k_true:.4f}")
    print(f"CLASS will get: H₀ = {H0_km_s_Mpc:.2f} km/s/Mpc  (from table at T_cmb)")
    print(f"CLASS will get: Ω_m (baryon) ≈ {omega_m_class:.6f}, Ω_r ≈ {omega_r_class:.6e}")
    print(f"Table: {n_loga} rows from a = {a_lock:.2e} to 1, T = {T_cmb_k} K at a=1")
    return {
        "eta": eta,
        "omega_k_true": omega_k_true,
        "a_lock": a_lock,
        "T_lock_gev": T_lock_gev,
        "H0_gev": H0_gev,
        "data": data,
        "outpath": outpath,
    }


if __name__ == "__main__":
    forward_4d_evolution()
