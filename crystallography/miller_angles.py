import tkinter as tk
from tkinter import filedialog
import numpy as np
import itertools

def extract_cell_and_spacegroup(cif_path):
    params = {}
    spacegroup_keys = [
        "_symmetry_space_group_name_H-M",
        "_space_group_name_H-M",
        "_symmetry_space_group_name_H-M_alt",
        "_space_group_name_H-M_alt"
    ]
    with open(cif_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("_cell_length_a"):
                params["a"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif line.startswith("_cell_length_b"):
                params["b"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif line.startswith("_cell_length_c"):
                params["c"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif line.startswith("_cell_angle_alpha"):
                params["alpha"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif line.startswith("_cell_angle_beta"):
                params["beta"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif line.startswith("_cell_angle_gamma"):
                params["gamma"] = float(line.split()[1].replace("(", "").replace(")", ""))
            elif any(line.startswith(key) for key in spacegroup_keys):
                value = line.split(maxsplit=1)[1].replace('"', '').replace("'", "")
                params["spacegroup"] = value
            elif line.startswith("_space_group_IT_number"):
                params["IT_number"] = line.split()[1]
    # Prompt for space group if not found
    if "spacegroup" not in params:
        print("Space group not found in CIF file.")
        params["spacegroup"] = input("Please enter the space group (H-M notation, e.g. 'P 21/c'): ").strip()
    # Check for other required parameters
    required = ["a", "b", "c", "alpha", "beta", "gamma", "spacegroup"]
    missing = [r for r in required if r not in params]
    if missing:
        print(f"Error: Missing required CIF parameters: {', '.join(missing)}")
        return None
    return params

def reciprocal_lattice(a, b, c, alpha_deg, beta_deg, gamma_deg):
    # Convert angles to radians
    alpha = np.deg2rad(alpha_deg)
    beta = np.deg2rad(beta_deg)
    gamma = np.deg2rad(gamma_deg)

    # Calculate volume of the unit cell
    volume = a * b * c * np.sqrt(
        1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2
        + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
    )

    # Direct lattice vectors in Cartesian coordinates
    a_vec = np.array([a, 0, 0])
    b_vec = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    c_vec = np.array([cx, cy, cz])

    # Reciprocal lattice vectors
    astar = 2 * np.pi * np.cross(b_vec, c_vec) / volume
    bstar = 2 * np.pi * np.cross(c_vec, a_vec) / volume
    cstar = 2 * np.pi * np.cross(a_vec, b_vec) / volume

    return astar, bstar, cstar

def plane_normal(hkl, astar, bstar, cstar):
    h, k, l = hkl
    return h * astar + k * bstar + l * cstar

def angle_between_planes(hkl1, hkl2, astar, bstar, cstar):
    n1 = plane_normal(hkl1, astar, bstar, cstar)
    n2 = plane_normal(hkl2, astar, bstar, cstar)
    dot = np.dot(n1, n2)
    norm1 = np.linalg.norm(n1)
    norm2 = np.linalg.norm(n2)
    if norm1 == 0 or norm2 == 0:
        return None
    cos_angle = np.clip(dot / (norm1 * norm2), -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))

def main():
    MAX_INDEX = 6
    MAX_RESULTS = 150
    # --- Step 1: Use Tkinter to select CIF file ---
    root = tk.Tk()
    root.withdraw()
    cif_path = filedialog.askopenfilename(
        title="Select a CIF file",
        filetypes=[("CIF files", "*.cif"), ("All files", "*.*")]
    )
    if not cif_path:
        print("No file selected. Exiting.")
        return

    print("Selected file path:", cif_path)

    # --- Step 2: Extract cell parameters and space group ---
    params = extract_cell_and_spacegroup(cif_path)
    if params is None:
        print("Please select a valid CIF file containing all required cell and space group information.")
        return

    print("\nExtracted parameters:")
    for k, v in params.items():
        print(f"{k}: {v}")

    # --- Step 3: Get user input for Miller plane and angle ---
    hkl_input = input("\nEnter Miller indices (h k l), e.g. 0 1 2: ").strip()
    hkl_ref = tuple(map(int, hkl_input.split()))
    target_angle = float(input("Enter target angle in degrees: ").strip())

    # --- Step 4: Calculate reciprocal lattice vectors ---
    astar, bstar, cstar = reciprocal_lattice(
        params["a"], params["b"], params["c"],
        params["alpha"], params["beta"], params["gamma"]
    )

    # --- Step 5: Compute angles to all planes up to (333) ---
    results = []
    for h, k, l in itertools.product(range(-MAX_INDEX, MAX_INDEX + 1), repeat=3):
        if (h, k, l) == (0, 0, 0):
            continue
        hkl = (h, k, l)
        angle = angle_between_planes(hkl_ref, hkl, astar, bstar, cstar)
        if angle is not None:
            results.append((hkl, angle))

    # --- Step 6: Sort and display ---
    results.sort(key=lambda x: abs(x[1] - target_angle))
    print(f"\nTop {MAX_RESULTS} planes sorted by proximity to {target_angle}° with respect to {hkl_ref}:")
    for hkl, angle in results[:MAX_RESULTS]:
        print(f"Plane {hkl}: {angle:.2f} degrees")

if __name__ == "__main__":
    main()