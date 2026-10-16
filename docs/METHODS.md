# Methods: what is implemented and what is assumed

This is a compact methods note intended to make the modeling choices understandable and reproducible. For the exact numerical implementation, please inspect the named Python modules.

## 1. Geometry and interpretation

All transport stages use a **one-dimensional idealized boundary layer**, not a full electrochemical cell. The cathodic side receives an imposed OH-equivalent flux, and the far boundary represents a bulk reservoir. The pH reported as 'near cathode' is generally at the **first finite-volume cell center**, not directly on the metal surface. The models do not include gas-bubble mixing, explicit fluid advection, detailed electrode kinetics, or membranes.

## 2. Reduced transport and pH (Stage 1)

The initial model uses a time-dependent diffusion equation for a reduced alkalinity variable:

\[
\frac{\partial A}{\partial t}=D\frac{\partial^2 A}{\partial x^2}.
\]

An imposed boundary flux represents cathodic alkalinity production. Carbonate buffering is represented by an instantaneous ideal-equilibrium calculation. This stage is a deliberately reduced model; it is verified numerically against an analytical diffusion case and through mesh refinement and mass balance. Files: `model.py`, `validation.py`, and `verification_run.py`.

## 3. Ionic migration (Stage 2)

For a species \(i\), the simplified Nernst–Planck flux combines diffusion and electric-field-driven migration:

\[
J_i=-D_i\frac{\partial c_i}{\partial x}
-\frac{z_i D_i F}{RT}c_i\frac{\partial\phi}{\partial x}.
\]

The calculation follows Na⁺, Cl⁻, H⁺ and OH⁻ with an electroneutrality constraint. This is a useful migration test in an idealized **unbuffered** electrolyte; it is not a model of an actual membrane-cell experiment. Files: `stage2_np.py`, `run_stage2.py`.

## 4. Carbonate equilibrium and scaling screen (Stage 3)

The subsequent model carries dissolved calcium, chloride, inorganic carbon, alkalinity and relevant ionic fluxes. Instantaneous water and carbonate equilibria are used to calculate local pH and idealized carbonate speciation. Sodium is constrained by charge balance in the model formulation.

The illustrative calcite saturation ratio is calculated from free-ion *concentrations* under an ideal approximation:

\[
\Omega_{\mathrm{ideal}}=\frac{[\mathrm{Ca}^{2+}][\mathrm{CO}_3^{2-}]}{K_{sp,\mathrm{ideal}}}.
\]

There is **no precipitation sink, solid phase, nucleation model or deposit mass**. An Ω value above one means supersaturation in the chosen thermodynamic approximation, not observed precipitation. Files: `stage3_carbonate.py`, `run_stage3.py`.

## 5. Activity effects and PHREEQC (Stage 5)

A separate approximate Davies calculation estimates nonideal activities. To check the solution chemistry more independently, PHREEQC v3 is run with the `phreeqc.dat` database using two different input decks:

- **Fixed-pH**: PHREEQC receives Stage 3 totals and the same imposed pH; this is the fairer comparison when isolating differences in speciation assumptions.
- **Charge-balanced pH**: PHREEQC is allowed to adjust pH to satisfy electrical neutrality; this also changes pH and therefore answers a different question.

PHREEQC reports \(SI_{\mathrm{calcite}}\), which is converted to saturation ratio using \(\Omega=10^{SI}\). No calcite equilibrium solid is imposed, since that would change dissolved inventories and make this specific comparison harder to interpret. Files: `stage5_thermo.py`, `phreeqc/benchmark_*.pqi`, `run_phreeqc.sh`, `compare_phreeqc.py`.

## 6. Spatial PHREEQC post-processing (Stage 6)

A Python script creates 160 PHREEQC input solutions from four modeled current conditions and 40 spatial **cell centers** per condition. PHREEQC then independently re-speciates the dissolved totals at imposed pH values. The post-processing plots ideal and PHREEQC Ω, total versus free calcium, dissolved CaCO₃(aq), charge-balance residuals and alkalinity differences.

This exchange is **one-way only**: the PHREEQC output does not change the ion fluxes computed in the transport stage. The resulting alkalinity and charge discrepancies are part of the current findings. Files: `stage6_spatial_prepare.py`, `stage6_spatial_compare.py`, `run_stage6.sh`.

## Verification versus validation

The Python test suite checks numerical and coding properties, including simplified analytic benchmarks, grid behavior, mass/charge accounting and output formatting. It does **not** verify that the modeled electrode chemistry describes a real experiment. The external PHREEQC calculation provides a useful thermodynamic *cross-check*, but is not independent experimental validation.

The PHREEQC input converts transport concentrations in mol/m³ (numerically equivalent to mmol/L) to approximately mmol/kg of water in dilute conditions; this is an approximation and should be checked for more concentrated solutions.
