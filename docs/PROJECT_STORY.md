# From a pH question to a more careful carbonate model

This document is the plain-English story behind the six stages. It is not a publication or a claim that the complete electrochemical system has been solved.

## Starting with a simple question

A cathode can produce hydroxide, making the liquid near its surface more alkaline than the rest of the cell. I wanted to see how far that pH change might extend, how it develops with time, and what happens when the liquid contains carbonate chemistry.

I began with a one-dimensional diffusion model, deliberately leaving out most of the complicated physics. It gave a useful first answer and, importantly, it could be checked against an analytical diffusion solution. That made it a good place to learn how numerical grid size and boundary conditions influence a result.

## Adding the physics that the first model missed

Diffusion is only one way ions move. In an electrolyte, electric fields can move positive and negative ions in different directions. The next step therefore used a simplified electroneutral Nernst–Planck description for Na⁺, Cl⁻, H⁺ and OH⁻. This stage helped separate concentration-driven transport from migration, while keeping the chemistry intentionally simple.

For a more realistic scaling question, the model then included carbonate equilibrium and calcium. At higher modeled pH, the ideal calculation returned much larger calcite saturation ratios. That was an interesting signal, but it was not evidence that solid calcite had formed: no nucleation, crystal growth, or solid-mass balance was included.

## Checking what other studies actually measured

Published studies of electrochemical desalination report pH shifts and carbonate-related scaling. But a pH measurement at a device outlet is not interchangeable with a modeled value within a few micrometres of an idealized cathode. Geometry, membranes, solution composition, flow, and electrode reactions all matter. I therefore treated the literature as context and a source of open questions—not as a direct validation dataset.

## Why the PHREEQC comparison mattered

Aqueous chemistry does not depend only on how many moles of a component are present. Activities and dissolved ion pairs matter, especially for doubly charged species such as calcium and carbonate.

The key comparison in the project kept the modeled pH fixed and asked PHREEQC to re-speciate the same component totals. Under the illustrative 0.5 A/m² case, the ideal indicator was Ω≈155.92; PHREEQC at the same pH gave Ω≈39.99. PHREEQC also assigned a substantial share of calcium to the dissolved neutral complex CaCO₃(aq).

I separately checked the *charge-balanced* PHREEQC calculation, where pH was allowed to change. That is a different question, so its result should not be described as only an activity correction.

## Exploring the whole boundary layer

The final stage wrote PHREEQC inputs for 40 grid-cell centers and four imposed-current scenarios. That made it possible to inspect how the chemical interpretation changes with distance from the cathode, rather than at one point alone. The user ran the external program on Ubuntu and reviewed the generated spatial plots.

This step also revealed a limitation worth keeping in view: PHREEQC and the Stage 3 transport model do not use fully consistent activity and alkalinity conventions. The one-way analysis exposes charge and alkalinity mismatches. It does not remove them.

## What I take away

The most useful outcome was not a single large saturation number. It was learning how a plausible-looking result can shift when the chemical assumptions are changed, and how numerical verification is different from experimental validation.

A rigorous next iteration would reconcile conserved component totals and alkalinity definitions, feed nonideal speciation back into the transport model, and test the predictions against local-pH measurements taken under comparable conditions. Only then would quantitative claims about real systems be appropriate.
