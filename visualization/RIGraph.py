import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D

# --- Data Setup ---
data = {
    'Wavelength (nm)': ['400','443.4','446','450','492','500','519.4','538','550','584','600','630','650','652.9'],
    'SiC': [2.759,0,0,2.716,0,2.684,0,0,2.66,0,2.643,0,2.629,0],
    'LiNbO$_3$': [2.439,0,2.382,0,2.345,0,0,2.32,0,2.302,0,2.287,0,2.276],
    '0428': [0,2.2231,0,0,0,0,2.1294,0,0,0,0,0,0,2.0587],
    # '0394': [0,??,0,0,0,0,2.15,0,0,0,0,0,0,??],
    # '0384': [0,??,0,0,0,0,2.23,0,0,0,0,0,0,??],
    '0367': [0,2.69,0,0,0,0,2.38,0,0,0,0,0,0,2.15],
    '0366': [0,2.71,0,0,0,0,2.43,0,0,0,0,0,0,2.25],
    # '0363': [0,??,0,0,0,0,2.27,0,0,0,0,0,0,??]
    '0356': [0,2.65,0,0,0,0,2.42,0,0,0,0,0,0,2.29],
    '0350': [0,2.66,0,0,0,0,2.39,0,0,0,0,0,0,2.23],
    '0348': [0,2.75,0,0,0,0,2.51,0,0,0,0,0,0,2.33],
    '0344': [0,2.17,0,0,0,0,1.82,0,0,0,0,0,0,1.58],
    '0334': [0,2.09,0,0,0,0,1.68,0,0,0,0,0,0,1.57],
    '0259': [0,2.85,0,0,0,0,2.64,0,0,0,0,0,0,2.38],
    # '0084': [0,??,0,0,0,0,2.68,0,0,0,0,0,0,2.34],
    '0029': [0,2.7,0,0,0,0,2.49,0,0,0,0,0,0,2.29],
    '0026': [0,2.55,0,0,0,0,2.37,0,0,0,0,0,0,2.26],
    '0011': [0,2.46,0,0,0,0,2.30,0,0,0,0,0,0,2.20],
    '2-ON': [0,2.18,0,0,0,0,2.09,0,0,0,0,0,0,2.05],
    'Pyrene': [0,2.38,0,0,0,0,2.21,0,0,0,0,0,0,2.13],
    '3-PB': [0,1.92,0,0,0,0,1.88,0,0,0,0,0,0,1.85]
}
df = pd.DataFrame(data)
df.iloc[:, 1:] = df.iloc[:, 1:].replace(0, np.nan)
df.set_index('Wavelength (nm)', inplace=True)

annotation_cfg = {
    'SiC': {'enabled': True, 'position': 'left'},
    'LiNbO$_3$': {'enabled': True, 'position': 'left'},
    '0428': {'enabled': False},
    # '0394': {'enabled': False},
    # '0384': {'enabled': False},
    '0367': {'enabled': False},
    '0366': {'enabled': False},
    # '0363': {'enabled': False},
    '0356': {'enabled': False},
    '0350': {'enabled': False},
    '0348': {'enabled': False},
    '0344': {'enabled': False},
    '0334': {'enabled': False},
    '0259': {'enabled': True, 'position': 'right'},
    # '0084': {'enabled': True, 'position': 'right'},
    '0029': {'enabled': True, 'position': 'right'},
    '0026': {'enabled': True, 'position': 'right'},
    '0011': {'enabled': True, 'position': 'right'},
    '2-ON': {'enabled': True, 'position': 'right'},
    'Pyrene': {'enabled': True, 'position': 'right'},
    '3-PB': {'enabled': True, 'position': 'right'},
}

color_map = {
    'SiC': "#DAB61A",
    'LiNbO$_3$': '#DAB61A',
    '0428': '#0064E0',   # Meta Blue
    # '0394': '#0064E0',
    # '0384': '#0064E0',
    '0367': '#0064E0',
    '0366': '#0064E0',
    # '0363': '#0064E0',
    '0356': '#0064E0',
    '0350': '#0064E0',
    '0348': '#0064E0',
    '0344': '#0064E0',
    '0334': '#0064E0',
    '0259': "#791717",   # Maroon
    # '0084': '#791717',
    '0029': '#791717',
    '0026': '#791717',
    '0011': '#791717',
    '2-ON': "#007A73",
    'Pyrene': '#791717',
    '3-PB': '#791717',
}
color_labels = {
    '#DAB61A': 'Inorganic',
    '#0064E0': 'Meta',
    '#791717': 'CCDC',
    '#007A73': 'Merck'
}

# --- Cauchy Fit Function ---
def cauchy_eq(lam, A, B, C):
    return A + B / lam**2 + C / lam**4

wavelengths = np.array(df.index, dtype=float)

plt.figure(figsize=(10,6))

# --- Plot Data and Cauchy Fit ---
for column in df.columns:
    valid_data = df[column].dropna()
    x = np.array(valid_data.index, dtype=float)
    y = np.array(valid_data.values, dtype=float)
    color = color_map.get(column, 'black')
    plt.scatter(x, y, color=color, s=0, label=f'{column} data' if column not in annotation_cfg or not annotation_cfg[column].get('enabled', False) else "")
    
    # Annotations
    cfg = annotation_cfg.get(column, {})
    if cfg.get('enabled', False) and len(x) > 0:
        pos = cfg.get('position', 'right')
        if pos == 'right':
            xy_idx = -1
            ha = 'left'
            xoffset = 10
        else:
            xy_idx = 0
            ha = 'right'
            xoffset = -10
        plt.annotate(
            column,
            xy=(x[xy_idx], y[xy_idx]),
            xytext=(xoffset, 0),
            textcoords='offset points',
            ha=ha, va='center',
            fontsize=12, fontweight='normal'
        )
    # --- Cauchy Fit Overlay ---
    if len(x) >= 3:
        try:
            p0 = [2.5, 1e4, 1e7]
            popt, _ = curve_fit(cauchy_eq, x, y, p0=p0)
            lam_fit = np.linspace(min(x), max(x), 200)
            y_fit = cauchy_eq(lam_fit, *popt)
            # Use dashed line for SiC and LiNbO$_3$, solid for others
            if column in ['SiC', 'LiNbO$_3$']:
                linestyle = '--'
            else:
                linestyle = '-'
            plt.plot(lam_fit, y_fit, linestyle=linestyle, color=color, alpha=1, linewidth=2)
            print(f"{column} Cauchy fit: A={popt[0]:.4f}, B={popt[1]:.2e}, C={popt[2]:.2e}")
        except Exception as e:
            print(f"Could not fit {column}: {e}")

# --- Custom Legend Handles ---
legend_handles = [
    Line2D([0], [0], color=color, lw=2, label=label)
    for color, label in color_labels.items()
]

plt.legend(handles=legend_handles, loc='lower left', frameon=False)
plt.text(0.9, 0.01, '*Values may increase due to measurement of non-optimised orientational axes.', 
         transform=plt.gca().transAxes, ha='right', fontsize=8)
plt.title('OSC RI dipersion curves (Cauchy fit)', fontsize=14, fontweight='normal')
plt.xlabel('Wavelength (nm)', fontsize=14, fontweight='normal')
plt.ylabel('RI - measured', fontsize=14, fontweight='normal')
plt.xlim(355, 690)
plt.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.7)
plt.tight_layout()
plt.show()
