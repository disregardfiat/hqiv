# HQIV Sandbox

Background integration for HQIV cosmology. Run `hqiv_background.py` to get age, proper time at z=14, d_A(z=14), and `hqiv_Ha.txt` for CLASS.

## Setup (recommended)

Use a **fresh venv with scipy** so the ODE integrator is stable:

```bash
# From repo root or anywhere
python -m venv hqiv_env
source hqiv_env/bin/activate   # Windows: hqiv_env\Scripts\activate
pip install numpy scipy matplotlib classy
cd sandbox
python hqiv_background.py
```

Takes ~5–10 min to install on decent hardware. If you get SSL errors on pip, try another network or install numpy/scipy via your OS (e.g. `apt install python3-scipy`).

## If numbers look wrong

- **Age / time / d_A way too large or tiny**  
  The horizon term and ODE are sensitive to units and the value of β. With **only numpy** (no scipy), the script uses a simple RK4 that can be unstable; install **scipy** and re-run so `solve_ivp` is used.
- **CLASS / full pipeline**  
  See main [README](../README.md) and `class_integration/MODIFICATIONS.md` for CLASS from source and HiCLASS-style custom background.

## Output

- Console: universe age (Gyr), proper time at z=14 (Myr), d_A(z=14) (Gpc), H0.
- `hqiv_Ha.txt`: columns `a`, `H_over_H0` for use as a custom background table in CLASS.
