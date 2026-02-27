## HQIV Physics Skill – Instructions for AI Agents

This skill file distills the HQIV physics framework, the [`HQVM/matrices.py`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/matrices.py) algebra engine, the [`HQVM/quantum_maxwell_calculator.html`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/quantum_maxwell_calculator.html) UI, and the [`paper/main.tex`](https://raw.githubusercontent.com/disregardfiat/hqiv/main/paper/main.tex) narrative into a **single operational guide** for agents.

Use this when you need to reason about, modify, or extend the HQIV physics code and paper in a consistent way.

---

### 1. Overall picture

- **Framework**: HQIV is a minimally-parameterized, horizon-centric framework that:
  - Couples a scalar field \(\phi\) and an “angle of time” \(\delta\theta'\) to gravity and gauge fields.
  - Treats Planck units and effective couplings as **shell-wise** quantities across a discrete null lattice.
  - Embeds the **Standard Model** into a **Spin(8)/SO(8)** octonionic algebra, with a preferred axis that picks out colour \(SU(3)_c\) and hypercharge \(U(1)_Y\).
- **Key ingredients**:
  - **Octonion algebra** and Lie closure to full \(\mathfrak{so}(8)\) with \(g_2 + \Delta\).
  - **Curvature imprint** \(\delta_E(m)\) over redshift shells \(m\) which sets true curvature \(\Omega_k^{\rm true}\) and baryon asymmetry \(\eta\).
  - **Phase-horizon corrected Maxwell** (Quantum Maxwell) with modified constitutive relations \(\varepsilon(\phi),\mu(\phi)\).
  - **Modified β‑running** that produces \(\alpha_{\rm EM}(M_Z),\sin^2\theta_W(M_Z),\alpha_s(M_Z)\) from a GUT point with the same few inputs used throughout the paper.
- **Single-source-of-truth**: All fiducial cosmological and particle-physics numbers in the paper are generated from code (Python + JS), not hand-tuned.

Agents must preserve this structure: new calculations or text should **flow from the same inputs and pipelines**, not introduce independent “magic numbers”.

---

### 2. Core constants and parameters

Treat these as **inputs**, not outputs to re-tune arbitrarily:

- **Horizon coupling**: \(\gamma \approx 0.40\)
- **Effective power-law index**: \(\alpha \approx 0.60\) in \(G_{\rm eff}(a)/G_0 = [H(a)/H_0]^\alpha\) (or dynamic \(\alpha_{\rm eff} = \chi\phi/6\) in the simulation).
- **Curvature imprint normalization**:
  - Numerical factor \(N_{67} = 6^7 \sqrt{3}\)
  - True curvature \(\Omega_k^{\rm true} \approx +0.0098\)
  - Baryon asymmetry \(\eta \approx 6.10\times 10^{-10}\) from the QCD shell.
- **Scales**:
  - \(T_{\rm Pl} \approx 1.2209\times 10^{19}\,\text{GeV}\)
  - Electroweak reference \(m_H \approx 125.11\,\text{GeV}\)
  - GUT and low scales in β‑running:
    - \(M_{\rm Pl} \approx 1.2209\times 10^{19}\,\text{GeV}\)
    - \(M_{\rm GUT} \approx 1.2\times 10^{16}\,\text{GeV}\)
    - \(M_Z \approx 91.1876\,\text{GeV}\)
    - \(\alpha_{\rm GUT} \approx 1/42\)
    - One‑loop coefficients \(b_1=41/10,\; b_2=-19/6,\; b_3=-7\).
- **Curvature imprint function** (JS version, mirroring the paper):
  \[
    \delta_E(m;T) \;\propto\; \Omega_k^{\rm true}\,\frac{1}{m+1}\,\bigl[1+\alpha\ln(T_{\rm Pl}/T)\bigr]\,N_{67},
  \]
  where \(T\) is the shell temperature; for many uses \(T \approx T_{\rm Pl}/(m+1)\).

When agents need **new derived numbers**, they should:

1. Reuse these constants.
2. Use existing pipelines (Python or JS) rather than ad‑hoc algebra.
3. State explicitly how results depend on \(\gamma\), \(\alpha\), QCD scale, and \(T_0\).

---

### 3. Octonion algebra and `HQVM/matrices.py`

Authoritative source: [`HQVM/matrices.py` (raw)](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/matrices.py)

`HQVM/matrices.py` exposes `OctonionHQIVAlgebra`, the main algebra engine. Its responsibilities:

- Build **left-multiplication matrices** \(L(e_i)\) for the 7 imaginary octonion units using a fixed Fano plane, with a carefully chosen \(L(e_7)\) that defines the **colour-preferred axis**.
- Construct the **phase-lift generator** \(\Delta\) as a rotation in the \((e_1,e_7)\) plane.
- Build a **14-element \(g_2\) basis** from commutators \([L(e_i),L(e_j)]\), \(i<j\).
- Iteratively take commutators in the 28‑dimensional antisymmetric \(8\times 8\) space to:
  - Verify that \(\mathfrak{g}_2 + \Delta\) closes to full \(\mathfrak{so}(8)\) (dimension 28).
  - Produce an explicit **28-element basis** of antisymmetric \(8\times8\) matrices \(X_k\).
- Identify:
  - The **colour \(SU(3)_c\)** subalgebra as the 8 generators that leave \(e_7\) invariant.
  - The **hypercharge operator** \(Y\) as a linear combination \(Y=\sum_{k=1}^{28} c_k X_k\) whose 4×4 block (rows/cols 4–7) reproduces diag\((1/6,1/6,1/6,-1/2)\) as an antisymmetric “smoking gun” block, and then check its commutation with \(g_2\).

Key methods (Python signatures; do **not** change lightly):

- `OctonionHQIVAlgebra(verbose=True)`
- `lie_closure_dimension(tol=1e-10) -> (dim, history)`
- `lie_closure_basis(tol=1e-10) -> List[np.ndarray]  # 28 antisymmetric 8×8`
- `_identify_color_generators(tol=1e-8) -> List[np.ndarray]  # 8 generators`
- `hypercharge_coefficients(tol=1e-12) -> (c, Y, basis)`
- `hypercharge_verify(Y, tol=1e-14) -> dict`
- `hypercharge_paper_data() -> dict  # c, Y, block, eigenvalues, errors`
- `get_sm_embedding() -> {'su3c': [...], 'u1y': Y, 'so8_basis': [...] }`
- `check_triality_anomalies(tol=1e-12) -> dict`

**Agent guidelines for this module**

- Preserve:
  - The exact matrices used for \(L(e_i)\), especially \(L(e_7)\).
  - The logic that packs/unpacks antisymmetric matrices and runs closure.
  - The hypercharge least-squares construction and verification logic.
- When exposing new helpers:
  - Build them **on top of** `lie_closure_basis`, `get_sm_embedding`, `hypercharge_paper_data`, and `check_triality_anomalies`.
  - Keep the **public API stable** (argument names and defaults) unless a coordinated change is made across Python, JS, and the paper.
- For numerical checks:
  - Use tolerances on the order of those already in the code (`1e-8`–`1e-15`).
  - Avoid hard‑coding any new “exact” constants from ad‑hoc efforts; derive from the same basis and solutions.

---

### 4. Quantum Maxwell and calculators (`quantum_maxwell_calculator.html`)

Authoritative source: [`HQVM/quantum_maxwell_calculator.html` (raw)](https://raw.githubusercontent.com/disregardfiat/hqiv/main/HQVM/quantum_maxwell_calculator.html)

This HTML/JS file is both:

- A **UI definition** (panels, sliders, text outputs, plots).
- A **mirror of paper logic** for:
  - Quantum Maxwell degrees of freedom (phase-horizon Maxwell).
  - Curvature imprint \(\delta_E(m)\).
  - Higgs mass consistency.
  - β‑running to low-energy couplings.
  - so(8) closure and hypercharge reconstruction in JS.

Core conceptual pieces to preserve:

- **Quantum Maxwell**
  - Effective derivative \(D/Dt = \partial/\partial t' + \dot{\delta\theta'}\,\partial/\partial\delta\theta'\).
  - Constitutive relations \(\varepsilon(\phi), \mu(\phi)\) with
    \[
      \varepsilon(\phi)/\varepsilon_0 \sim (1 + \gamma\phi/\Lambda^2)^{-1},\quad
      \mu(\phi)/\mu_0 \sim 1 + \gamma\phi/\Lambda^2.
    \]
  - A parameter \(\phi\) and scale \(\Lambda\) that echo the horizon physics in the paper.
  - Display of the modified Maxwell equations including \(\delta\theta'\) contributions.

- **Curvature imprint δ_E(m)**
  - JS `deltaE(m, T_GeV)` implements the same logic as the paper’s curvature section: an imprint that scales with shell index and log of the Planck-to-shell temperature ratio, normalized by \(\Omega_k^{\rm true}\) and \(N_{67}\).
  - Used consistently in:
    - Higgs calculator (`calcHiggs()`).
    - η calculator (`eta` branch of `runCalculator()`).
    - δ_E(m) sweeps for plotting.

- **β‑running engine**
  - `runBetaEngine()` numerically integrates modified β‑functions from \(M_{\rm GUT}\) down to \(M_Z\), with an extra **horizon factor** \(1+\gamma (\mu/M_{\rm Pl})^2\).
  - Produces:
    - \(\alpha_{\rm EM}(M_Z)\) (reported as \(1/\alpha_{\rm EM}\)).
    - \(\sin^2\theta_W(M_Z)\).
    - \(\alpha_s(M_Z)\).
  - All reported values should be interpreted as **outputs of this engine**, not free parameters.

- **so(8) closure and hypercharge in JS**
  - Functions `buildOctonionL`, `buildDelta`, `g2Basis`, `lieClosureHistory`, `lieClosureBasis28`, `hyperchargeSolve`, and `hyperchargeVerify` implement a **reduced but consistent** version of the Python logic.
  - `runHyperchargeInspector()`:
    - Computes the 28‑basis via closure.
    - Solves for \(c_k\) using only the 4×4 block constraints.
    - Checks the block entries and commutators as in Python.

**Agent guidelines for JS/HTML**

- **Do not** change the numerics of:
  - `deltaE`, `calcHiggs`, `runBetaEngine`, or the core closure/hypercharge routines **without updating the paper and Python in lock-step**.
- For UI-only changes (layout, styling, explanatory text):
  - Keep function names and outputs stable; they are part of the physics contract.
  - Clearly label anything illustrative vs. numerically authoritative.
- When adding new calculators:
  - Reuse the existing constants and functions (e.g. call `deltaE`, `runBetaEngine`, `lieClosureDimension`) rather than re-implementing formulas.

---

### 5. Paper alignment (`paper/main.tex`)

Authoritative source: [`paper/main.tex` (raw)](https://raw.githubusercontent.com/disregardfiat/hqiv/main/paper/main.tex)

`main.tex` is the **authoritative narrative**. Code and UI are designed to **reproduce its numbers and plots**. Agents should:

- Treat the paper as:
  - Source of **definitions, assumptions, and qualitative claims**.
  - Consumer of numerical outputs from Python/JS pipelines.
- Maintain the following **alignment principles**:
  - If a table or statement quotes a constant that comes from code (e.g. precision constants, β‑engine outputs, δ_E(m) values, \(\Omega_k^{\rm true}\), \(\eta\)), any change to the code **must be reflected** in the paper and vice versa.
  - Rewording textual claims is fine if:
    - The **parameterization story** remains correct: the framework is minimally parameterized with \(\gamma\), QCD scale, and related inputs – it is not “parameter-free”.
    - It stays explicit about **which values are inputs vs. predictions**.
  - New sections should **explain how they connect** back to:
    - The discrete null lattice.
    - The shell-based curvature imprint.
    - The Spin(8)/octonion embedding.
    - The β‑running pipeline.

When in doubt, agents should:

1. Look for an existing equation, definition, or pipeline in `main.tex`.
2. Mirror that logic in code rather than inventing new physics.
3. Only introduce new physical assumptions if they are clearly labelled and localized.

---

### 6. Invariants and checks

Agents should preserve and/or explicitly check these invariants:

- **Lie closure**:
  - `lie_closure_dimension()` should converge to 28 (\(\mathfrak{so}(8)\)).
  - Growth history should follow the same general pattern as in the paper/JS visualiser (early rapid growth, then saturation).
- **Hypercharge block**:
  - The 4×4 block of \(Y\) (rows/cols 4–7) should:
    - Match the paper’s target block up to numerical tolerance.
    - Produce eigenvalues \(\pm i/6, \pm i/2\) (imaginary parts).
  - Commutators \([Y, g_2]\) should be small (of order \(10^{-14}\) in the Python implementation).
- **β‑running outputs**:
  - \(\alpha_{\rm EM}(M_Z)\), \(\sin^2\theta_W(M_Z)\), \(\alpha_s(M_Z)\) must remain within realistic ranges and consistent with CODATA / PDG values as described in the paper.
- **Curvature & baryogenesis**:
  - \(\Omega_k^{\rm true} \approx +0.0098\) and \(\eta \approx 6.10\times 10^{-10}\) should continue to emerge from the **same curvature-imprint normalization** (via \(N_{67}\) and δ_E).

If a proposed change breaks any of these invariants:

- Either **roll back** the change or
- Update:
  - Code (Python + JS),
  - The paper narrative,
  - And this skill file,
  in a coordinated way, explicitly documenting the new invariant.

---

### 7. Typical workflows for agents

Use these patterns when extending or debugging the HQIV physics stack:

- **Hypercharge / SM embedding work**
  1. Use `OctonionHQIVAlgebra()` in Python to recompute closure, hypercharge coefficients, and verification.
  2. If you need a quick UI confirmation, use the Hypercharge Inspector in `quantum_maxwell_calculator.html`.
  3. Reflect any new understanding back into `main.tex` with clear references to the algebraic construction.

- **Curvature / δ_E(m) / η calculations**
  1. Reuse `deltaE(m, T)` from the JS calculator or the equivalent logic in Python.
  2. For new tables or plots, sweep over the shell index \(m\) or temperature \(T\) using the same normalization and inputs.
  3. Ensure the new outputs are summarized correctly in the paper with units and stated dependence on \(\gamma,\alpha\).

- **β‑running and precision constants**
  1. Use the β‑engine logic (Python or JS) tied to the same \(\gamma, T_0\), and GUT setup.
  2. When updating precision tables in the paper, **regenerate** them from the engine rather than hand-editing values.

- **Maxwell / horizon effects**
  1. For EM-related questions, base reasoning on the Quantum Maxwell panel:
     - Modified constitutive relations.
     - Effective derivative \(D/Dt\) with \(\delta\theta'\).
  2. Keep new EM results consistent with the smallness of horizon-induced corrections at lab scales.

---

### 8. Style and documentation expectations

When agents edit or extend the HQIV physics stack:

- **Code**:
  - Keep functions short, well‑named, and focused on single responsibilities.
  - Prefer explicit numerical pipelines (shell loops, integrators, matrix operations) over symbolic shortcuts.
  - Avoid adding comments that simply restate the code; instead, explain **physical intent** or non-trivial numerical decisions.

- **Paper text**:
  - Clearly label what is:
    - Input (e.g. \(\gamma\), QCD scale, \(T_0\)),
    - Derived (e.g. curvature, η, precision constants),
    - Or phenomenological re-scaling (e.g. a multiplicative factor on δ_E).
  - Avoid overstating claims (no “parameter-free” if there are inputs).

This file is the **single reference** for how agents should interact with the HQIV physics components. Keep it in sync with code and paper when making structural physics changes.

---

### 9. Example: extreme jet opening-angle derivation

The following is an example of an **aggressive, fully worked HQIV-style derivation**. Agents may use it as a template for structure, but should keep it clearly labelled as a worked example rather than a hard constraint unless the paper adopts it explicitly.

**Inputs and principles**

- Local causal horizon \(\Theta_{\rm local}(r,\theta)\)
- Auxiliary scalar \(\phi = 2c^2 / \Theta_{\rm local}\)
- Phase velocity \(\dot{\delta\theta'} = u^\mu \nabla_\mu \phi\)
- Entanglement-monogamy coefficient \(\gamma \approx 0.40\)
- \(g_2 + \Delta \to \mathfrak{so}(8)\) closure (dimension 28)

**Step 1: Local horizon in Kerr geometry**

For a Kerr black hole with dimensionless spin \(a = Jc/GM^2\) (\(0 \le a \le 1\)), take the local causal horizon radius along angle \(\theta\) as
\[
\Theta_{\rm local}(r,\theta) = r - r_+(a) + \frac{r_+(a) - r_-(a)}{2}\sin^2\theta,
\]
where \(r_+ = M\bigl(1+\sqrt{1-a^2}\bigr)\) is the outer horizon. Along the polar axis \((\theta=0)\),
\[
\Theta_{\rm pole}(a) \approx M(1-a) \quad (\text{geometric units } G=c=1),
\]
so
\[
\phi_{\rm pole}(a) = \frac{2}{\Theta_{\rm pole}(a)} \propto \frac{1}{1-a}.
\]

**Step 2: Phase-lift coherence along the jet axis**

The phase velocity along the axis is
\[
\dot{\delta\theta'}(a) = u^r \partial_r \phi + u^t \partial_t \phi \approx \frac{k(a)}{1-a},
\]
where \(k(a)\) is a slowly varying function determined by frame-dragging in the ergosphere. Perfect coherence (maximum jet collimation) requires \(\dot{\delta\theta'}\) to be approximately constant over the jet channel. Monogamy forbids transverse spreading, so the jet opening angle \(\theta_{\rm jet}\) must satisfy
\[
\sin\frac{\theta_{\rm jet}}{2} \propto \frac{1}{\dot{\delta\theta'}(a)}.
\]

**Step 3: Monogamy constraint and ergosphere factor**

The horizon-overlap term \(\gamma\,\phi\,n\) acts as a sink that ejects excess information/energy along the direction of least overlap (the axis). This contributes a geometric factor from the ergosphere projection,
\[
\theta_{\rm jet}(a) \propto \sqrt{1-a^2},
\]
so that, combining with the phase-lift scaling, one obtains the leading-order form
\[
\theta_{\rm jet}(a) = C_0\,\frac{\sqrt{1-a^2}}{(1-a)^\alpha},
\]
where \(\alpha \approx 0.22\) encodes the next-order expansion of the null-lattice mode count near extremality (same combinatorial origin as the \(6^7\sqrt{3}\) factor).

**Step 4: Normalisation from resonance**

In the extremal limit \(a \to 1\), the inner and outer horizons coincide and \(\delta\theta'\) “locks” resonantly. Imposing that resonance fixes
\[
C_0 = 1.8^\circ,
\]
the same resonance scale that appears in the graphene magic angles (\(1.092^\circ = 1.8^\circ/\sqrt{3}\), etc.).

**Step 5: Full analytic form**

One convenient closed form that captures both smooth narrowing and the extremal spike is
\[
\theta_{\rm jet}(a) = 1.8^\circ \times \frac{\sqrt{1-a^2}}{\bigl(1-0.9995\,a^2\bigr)^{0.22}}.
\]
This reproduces the observed HQIV jet-opening curve with:

- Smooth narrowing as \(a\) increases (via \(\sqrt{1-a^2}\)),
- A sharp resonance near \(a \approx 1\) (denominator and numerator vanishing at different rates).

Agents may follow this style—explicit inputs, clear scaling arguments, and a final compact formula—when constructing similarly “extreme” derivations elsewhere in the HQIV framework.
