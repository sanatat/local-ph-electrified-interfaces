# Contribution, scientific scope and open questions

## What this project contributes

This project combines a series of **small, testable models** with an independent equilibrium-chemistry check. The useful contribution is the *workflow*: moving from near-cathode diffusion to ion migration and carbonate equilibria, then checking how a more complete speciation model changes the interpretation.

For the illustrative 0.5 A/m² case, the idealized carbonate calculation returned Ω≈155.92, while PHREEQC at the **same pH** returned Ω≈39.99. That contrast provided a concrete reason to study activity corrections and calcium complexation. A second useful outcome was finding nonzero charge and alkalinity differences in the one-way coupling instead of quietly treating the combined workflow as self-consistent.

## Is it novel?

**A claim of scientific novelty has not been established.** Nernst–Planck transport, carbonate equilibrium, PHREEQC and their combinations all have substantial prior literature. This repository does not claim to have discovered a new chemical mechanism or invented a previously unknown algorithm.

What may distinguish this as a **research portfolio** is the clear progression of models, quantitative checks, fixed-pH versus charge-balanced comparisons, and transparent documentation of mismatches. That is an educational and methodological strength, not proof of priority over existing publications.

## What the current results do not show

| What is present | What it should not be mistaken for |
|---|---|
| 33 successful coding and numerical checks | Laboratory validation of the pH profile |
| An idealized first-cell pH calculation | A measurement at the electrode surface |
| Ω greater than one | A measured mass or rate of solid CaCO₃ deposition |
| CaCO₃(aq) formation | Precipitated CaCO₃(s) |
| PHREEQC calculation at imposed pH | A closed-loop activity-corrected transport solution |
| 160 spatial chemistry inputs | 160 independent experiments |

## The unresolved research gap

The next substantial task is to **align conservation variables and solution chemistry** between Stage 3 and PHREEQC. Both must use compatible definitions of total alkalinity, dissolved components, ionic activities, and charge balance. Only then would it be reasonable to create a coupled iterative transport/speciation solution.

Further additions—membrane charge and selectivity, flow, bubble formation, electrode kinetics and precipitation kinetics—should be driven by a specific experimental system. Adding all of them without matched data would increase model complexity without necessarily improving credibility.

## Proposed validation strategy

1. Establish agreement in a simple, explicitly defined equilibrium case with consistent pH and total-component inputs.
2. Diagnose and reduce charge/alkalinity residuals **without concealing them**.
3. Test grid and time-step dependence after chemistry is coupled back into transport.
4. Identify an independently measured local-pH dataset with comparable solution chemistry, geometry and applied current, and evaluate prediction error on held-out observations.
5. Only then explore deposition onset and scaling mitigation under appropriate precipitation kinetics.

## Relationship to published work

Studies by Arulrajan and colleagues investigate pH changes and carbonate scaling in electrochemical desalination (2021: https://doi.org/10.1021/acs.est.1c04479; 2024: https://doi.org/10.1016/j.jwpe.2024.105094). They provide a reason to ask this question, but their device-specific bulk/outlet measurements are **not interchangeable** with a computed first-cell cathode pH in this idealized system.

This repository is independent. References are scholarly background, not claims of supervision or collaboration.
