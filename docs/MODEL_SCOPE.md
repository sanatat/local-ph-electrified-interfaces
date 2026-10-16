# Model scope at a glance

This is an independent **computational exploration**, not an experimental study or a validated device model.

| Part | Included | Left out |
|---|---|---|
| Diffusion model | 1-D time-dependent transport and test cases | Electrode surface reactions and flow |
| Nernst–Planck model | Idealized Na⁺/Cl⁻/H⁺/OH⁻ diffusion and migration | Complex buffers and membrane physics |
| Carbonate model | Ideal carbonate equilibrium, calcium inventory, Ω screen | Nonideal complexes and solid precipitation |
| Davies + PHREEQC | Approximate activity checks and independent chemical speciation | Fully matched underlying thermodynamic equations |
| Spatial PHREEQC step | 160 fixed-pH input snapshots and diagnostic plots | Feedback from PHREEQC to the transport solver |

Remember: **Ω>1 is a thermodynamic tendency, not a deposition measurement**. The neutral species **CaCO₃(aq)** is dissolved calcium carbonate, not solid calcite. Results at the first cell center are not values measured on the metal surface.

The model uses dilute-solution unit conversions from mol/m³ to approximately mmol/kg water. Check this approximation before applying it to concentrated electrolytes.

For more context, read [Methods](METHODS.md), [Results](RESULTS.md), and [Limitations](SCIENTIFIC_SCOPE_AND_NOVELTY.md).
