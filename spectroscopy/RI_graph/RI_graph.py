import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import os

# --- Reference data (fixed, full wavelength range) ---
reference_data = {
    'SiC': {
        'wavelengths': [400, 450, 500, 550, 600, 650],
        'ri': [2.759, 2.716, 2.684, 2.66, 2.643, 2.629],
        'color': "#000000",
        'label_position': 'left'
    },
    'LiNbO$_3$': {
        'wavelengths': [400, 446, 492, 538, 584, 630],
        'ri': [2.439, 2.382, 2.345, 2.32, 2.302, 2.287],
        'color': "#000000",
        'label_position': 'left'
    }
}

# --- Load measured data from local CSV ---
script_dir = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(script_dir, 'data.csv'))

# --- Color mapping by Type ---
type_colors = {
    'Inorganic': "#000000",
    'Meta': '#0064E0',
    'CCDC': "#791764",
    'Merck': '#007A73',
}

wavelength_columns = ['443.4nm', '519.4nm', '652.9nm']
wavelengths = [443.4, 519.4, 652.9]

def cauchy_eq(lam, A, B, C):
    return A + B / lam**2 + C / lam**4

def is_valid_number(val):
    """Check if a value can be converted to a valid float."""
    try:
        f = float(val)
        return not np.isnan(f)
    except (ValueError, TypeError):
        return False

def get_valid_points(row):
    """Extract valid (wavelength, ri) pairs from a row."""
    valid_x = []
    valid_y = []
    for wl, col in zip(wavelengths, wavelength_columns):
        val = row[col]
        if is_valid_number(val):
            valid_x.append(wl)
            valid_y.append(float(val))
    return np.array(valid_x), np.array(valid_y)

plt.figure(figsize=(10, 6))

# --- Plot reference materials (dashed lines) ---
for material, data in reference_data.items():
    x = np.array(data['wavelengths'], dtype=float)
    y = np.array(data['ri'], dtype=float)
    color = data['color']
    
    try:
        p0 = [2.5, 1e4, 1e7]
        popt, _ = curve_fit(cauchy_eq, x, y, p0=p0)
        lam_fit = np.linspace(min(x), max(x), 200)
        y_fit = cauchy_eq(lam_fit, *popt)
        
        plt.plot(lam_fit, y_fit, linestyle='--', color=color, linewidth=2)
        print(f"{material} Cauchy fit: A={popt[0]:.4f}, B={popt[1]:.2e}, C={popt[2]:.2e}")
    except Exception as e:
        print(f"Could not fit {material}: {e}")
    
    pos = data['label_position']
    if pos == 'left':
        xy_idx, ha, xoffset = 0, 'right', -10
    else:
        xy_idx, ha, xoffset = -1, 'left', 10
    plt.annotate(material, xy=(x[xy_idx], y[xy_idx]), xytext=(xoffset, 0),
                 textcoords='offset points', ha=ha, va='center', fontsize=12)

# --- Plot measured materials ---
for _, row in df.iterrows():
    material = row['Material']
    
    # Skip if not included
    if row.get('Include', 'Yes') != 'Yes':
        print(f"Skipping {material}: Include = No")
        continue
    
    mat_type = row['Type']
    
    # Use color override if provided, otherwise use type color
    color_override = row.get('Color Override', '')
    if pd.notna(color_override) and str(color_override).strip():
        color = str(color_override).strip()
    else:
        color = type_colors.get(mat_type, 'black')
    
    # Get valid data points
    x, y = get_valid_points(row)
    num_points = len(x)
    
    if num_points == 0:
        print(f"Skipping {material}: No valid RI data")
        continue
    
    elif num_points == 1:
        # Single point - plot as marker with label
        plt.scatter(x, y, color=color, s=50, zorder=5)
        plt.annotate(material, xy=(x[0], y[0]), xytext=(10, 0),
                     textcoords='offset points', ha='left', va='center', fontsize=10)
        print(f"{material}: Single point at {x[0]}nm = {y[0]}")
    
    elif num_points == 2:
        # Two points - plot as straight line
        plt.plot(x, y, linestyle='-', color=color, linewidth=2, marker='o', markersize=4)
        print(f"{material}: Line between {x[0]}nm and {x[1]}nm")
        
        # Annotation
        if row.get('Show Label', 'No') == 'Yes':
            pos = row.get('Label Position', 'right')
            if pos == 'right':
                xy_idx, ha, xoffset = -1, 'left', 10
            else:
                xy_idx, ha, xoffset = 0, 'right', -10
            plt.annotate(material, xy=(x[xy_idx], y[xy_idx]), xytext=(xoffset, 0),
                         textcoords='offset points', ha=ha, va='center', fontsize=12)
    
    else:  # num_points == 3
        # Three points - Cauchy fit
        try:
            p0 = [2.5, 1e4, 1e7]
            popt, _ = curve_fit(cauchy_eq, x, y, p0=p0)
            lam_fit = np.linspace(min(x), max(x), 200)
            y_fit = cauchy_eq(lam_fit, *popt)
            
            plt.plot(lam_fit, y_fit, linestyle='-', color=color, linewidth=2)
            print(f"{material} Cauchy fit: A={popt[0]:.4f}, B={popt[1]:.2e}, C={popt[2]:.2e}")
        except Exception as e:
            print(f"Could not fit {material}: {e}")
            continue
        
        # Annotation
        if row.get('Show Label', 'No') == 'Yes':
            pos = row.get('Label Position', 'right')
            if pos == 'right':
                xy_idx, ha, xoffset = -1, 'left', 10
            else:
                xy_idx, ha, xoffset = 0, 'right', -10
            plt.annotate(material, xy=(x[xy_idx], y[xy_idx]), xytext=(xoffset, 0),
                         textcoords='offset points', ha=ha, va='center', fontsize=12)

# --- Legend ---
legend_handles = [Line2D([0], [0], color=c, lw=2, label=l) for c, l in 
                  [('#DAB61A', 'Inorganic'), ('#0064E0', 'Meta'), 
                   ('#791717', 'CCDC'), ('#007A73', 'Merck')]]

plt.legend(handles=legend_handles, loc='lower left', frameon=False)
plt.text(0.9, 0.01, '*Values may increase due to measurement of non-optimised orientational axes.', 
         transform=plt.gca().transAxes, ha='right', fontsize=8)
plt.title('OSC RI dispersion curves (Cauchy fit)', fontsize=14)
plt.xlabel('Wavelength (nm)', fontsize=14)
plt.ylabel('RI - measured', fontsize=14)
plt.xlim(355, 690)
plt.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.7)
plt.tight_layout()
plt.show()