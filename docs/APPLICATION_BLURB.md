# Describing this work in a research application

These drafts are designed to sound like a real researcher reflecting on a specific problem. **Use them only if you have run, reviewed and can explain the models.** Adjust the wording to reflect your actual contribution; do not claim that you independently wrote every line of code if the code was prepared with assistance.

## Short CV description

**Independent computational study — Local pH and carbonate speciation (2026).** Reproduced and examined a staged Python finite-volume model for cathode-adjacent pH and ion migration, and used PHREEQC to compare ideal and activity-aware calcite saturation estimates. Worked with a 160-case spatial speciation workflow and documented charge-balance and alkalinity limitations. Numerical study; no laboratory validation.

## Motivation-letter paragraph

I became interested in the way electrical current changes chemistry close to an electrode, where the conditions can be quite different from those in the bulk solution. To explore this, I worked with a reproducible Python and PHREEQC model of ionic transport, carbonate buffering and calcite saturation. What stood out to me was how much the saturation estimate changed when I moved from an ideal-concentration calculation to a more complete equilibrium model. The remaining charge-balance and alkalinity differences also showed me why careful experimental validation matters. I would like to build on this foundation by combining more rigorous transport modeling with local-pH measurements in a research group working on electrically controlled chemistry.

## GitHub description (short)

`Python + PHREEQC exploration of local cathodic pH, ionic transport, carbonate chemistry and calcite saturation; 33 numerical tests. Not experimentally validated.`

## Suggested repository topics

`electrochemistry`, `nernst-planck`, `phreeqc`, `carbonate-chemistry`, `reactive-transport`, `scientific-python`, `reproducible-research`
