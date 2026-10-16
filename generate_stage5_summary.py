"""Creates explanatory experimental figure from article-text numbers. It is NOT model fit."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
x=pd.read_csv(p/'experimental_reference_points.csv')
fig,ax=plt.subplots(figsize=(8.4,4.7))
ax.plot(x['current_magnitude_A_m2'],x['reported_ph_max'],'o-',lw=2)
for xx,yy in zip(x['current_magnitude_A_m2'],x['reported_ph_max']):ax.annotate(f'{yy:.1f}',(xx,yy),xytext=(0,8),textcoords='offset points',ha='center')
ax.set_ylim(5.7,7.8)
ax.set_title('Published FLIM pH maxima: flowing CO₂ electrolyzer (not model validation)')
ax.set_xlabel('Magnitude of cathodic current density (A/m²)')
ax.set_ylabel('Reported catholyte pHmax (not electrode-surface pH)')
ax.grid(alpha=.3)
fig.text(.5,.015,'Baumgartner et al. (2023), DOI:10.1021/acssuschemeng.3c01773. N₂-purged, Re=5.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(p/'results_stage5'/'published_flm_unmatched.png',dpi=180);plt.close(fig)
