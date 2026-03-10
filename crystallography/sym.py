import os
import tkinter as tk
from tkinter import filedialog
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter

def process_cifs_in_folder():
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(title="Select folder containing CIF files")
    if not folder_path:
        print("No folder selected, exiting.")
        return

    output_folder = os.path.join(folder_path, "sym_fixed")
    os.makedirs(output_folder, exist_ok=True)

    symprec = 0.5  # Symmetry tolerance

    for filename in os.listdir(folder_path):
        if filename.lower().endswith(".cif"):
            input_path = os.path.join(folder_path, filename)
            print(f"Processing {filename}...")

            try:
                structure = Structure.from_file(input_path)
                sga = SpacegroupAnalyzer(structure, symprec=symprec)

                # Get refined structure with symmetry applied
                refined_structure = sga.get_refined_structure()

                output_path = os.path.join(output_folder, filename)
                # Write CIF with symmetry info and tolerance
                CifWriter(refined_structure, symprec=symprec).write_file(output_path)

                print(f"Written symmetrized CIF to {output_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    process_cifs_in_folder()
