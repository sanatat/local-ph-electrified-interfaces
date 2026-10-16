"""Prepare one-way fixed-pH PHREEQC speciation for 160 Stage-3 spatial snapshots."""
import csv
from dataclasses import replace
from pathlib import Path
from stage3_carbonate import CarbonateCase, solve_carbonate

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results_stage6'
PQI = ROOT / 'phreeqc' / 'spatial_fixed_ph.pqi'
CURRENTS = (0.0, 0.15, 0.5, 1.0)


def prepare():
    OUT.mkdir(exist_ok=True)
    PQI.parent.mkdir(exist_ok=True)
    snapshots = []
    base = CarbonateCase()
    index = 1
    for current in CURRENTS:
        case = replace(base, current_Am2=current)
        result = solve_carbonate(case, times=(0., case.time_s))
        for cell, distance in enumerate(result['x_um']):
            snapshots.append({
                'soln': index, 'current_A_m2': current,
                'time_s': case.time_s, 'distance_um': float(distance),
                'ph_stage3': float(result['pH'][-1, cell]),
                'na_total_mmolkg': float(result['na_molm3'][-1, cell]),
                'cl_total_mmolkg': float(result['cl_molm3'][-1, cell]),
                'ca_total_mmolkg': float(result['ca_molm3'][-1, cell]),
                'dic_total_mmolkg': float(result['dic_molm3'][-1, cell]),
                'ta_stage3_mmolkg': float(result['ta_molm3'][-1, cell]),
                'omega_ideal': float(result['omega_ideal'][-1, cell])})
            index += 1
    with (OUT / 'input_stage6_profiles.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(snapshots[0]))
        writer.writeheader()
        writer.writerows(snapshots)
    text = [
      'TITLE Stage 6: one-way spatial fixed-pH PHREEQC speciation',
      '# 160 cells (40 per current); transport inventories in mmol/m3 = mmol/L ≈ mmol/kg water',
      '# Keep model pH fixed, do NOT impose charge balancing or calcite EQUILIBRIUM_PHASES',
      '# CaCO3 (aqueous species) is not CaCO3(s) or precipitation.',
      'SELECTED_OUTPUT 1',
      '    -file results_stage6/phreeqc_selected.tsv',
      '    -reset false',
      '    -solution true',
      '    -alkalinity true',
      '    -percent_error true',
      '    -high_precision true',
      'USER_PUNCH 1',
      '    -headings ph_activity si_calcite ionic_strength ca_free_molkg co3_free_molkg caco3_aq_molkg cahco3_aq_molkg caoh_aq_molkg naco3_aq_molkg',
      '    -start',
      '    10 PUNCH -LA("H+"), SI("Calcite"), MU, MOL("Ca+2"), MOL("CO3-2"), MOL("CaCO3"), MOL("CaHCO3+"), MOL("CaOH+"), MOL("NaCO3-")',
      '    -end',
      'END']
    for s in snapshots:
        text.extend([f"# case current={s['current_A_m2']:g} A/m2 x={s['distance_um']:.3f} um t=180s",
           f"SOLUTION {s['soln']} case_{s['soln']}",
           '    temp 25', '    units mmol/kgw',
           f"    pH {s['ph_stage3']:.12g}", '    pe 4',
           f"    Na {s['na_total_mmolkg']:.12g}",
           f"    Cl {s['cl_total_mmolkg']:.12g}",
           f"    Ca {s['ca_total_mmolkg']:.12g}",
           f"    C(4) {s['dic_total_mmolkg']:.12g}",
           '    -water 1', 'END'])
    PQI.write_text('\n'.join(text) + '\n', encoding='utf-8')
    print(f'Prepared {len(snapshots)} spatial samples: {PQI} (fixed pH, no precipitation)')
    return snapshots

if __name__ == '__main__': prepare()
