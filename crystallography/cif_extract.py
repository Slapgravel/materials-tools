import os
import tkinter as tk
from tkinter import filedialog
import csv
import re

def parse_value_with_esd(value_str):
    """
    Parse a value that may contain an ESD in parentheses.
    e.g., "91.735(16)" -> ("91.735", "16")
         "90" -> ("90", "")
    """
    if value_str is None:
        return None, None
    
    match = re.match(r'^([\d.]+)\((\d+)\)$', value_str)
    if match:
        return match.group(1), match.group(2)
    else:
        return value_str, ""
    
def convert_esd(esd_val):
    if esd_val:
        return round(float(esd_val) * 0.0001, 4)
    return None

def extract_parameters(folder_path):
    parameters = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".cif"):
            filepath = os.path.join(folder_path, filename)
            with open(filepath, 'r') as file:
                lines = file.readlines()
                space_group_number = None
                a = None
                b = None
                c = None
                alpha = None
                beta = None
                gamma = None
                for line in lines:
                    if "_space_group_IT_number" in line:
                        space_group_number = line.strip().split()[-1]
                    elif "_cell_measurement_temperature" in line:
                        temp = line.strip().split()[-1]
                    elif "_cell_length_a" in line:
                        a = line.strip().split()[-1]
                    elif "_cell_length_b" in line:
                        b = line.strip().split()[-1]
                    elif "_cell_length_c" in line:
                        c = line.strip().split()[-1]
                    elif "_cell_angle_alpha" in line:
                        alpha = line.strip().split()[-1]
                    elif "_cell_angle_beta" in line:
                        beta = line.strip().split()[-1]
                    elif "_cell_angle_gamma" in line:
                        gamma = line.strip().split()[-1]
                    elif "_cell_volume" in line:
                        volume = line.strip().split()[-1]
                    elif "REM R1 =" in line:
                        Rone = line.strip().split()[3]
                        Rone = round(float(Rone)*100 , 2)
                
                # Parse each value to separate value and ESD
                temp_val, temp_esd = parse_value_with_esd(temp)
                a_val, a_esd = parse_value_with_esd(a)
                b_val, b_esd = parse_value_with_esd(b)
                c_val, c_esd = parse_value_with_esd(c)
                alpha_val, alpha_esd = parse_value_with_esd(alpha)
                beta_val, beta_esd = parse_value_with_esd(beta)
                gamma_val, gamma_esd = parse_value_with_esd(gamma)
                vol_val, vol_esd = parse_value_with_esd(volume)
                
                a_esd = convert_esd(a_esd)
                b_esd = convert_esd(b_esd)
                c_esd = convert_esd(c_esd)
                alpha_esd = convert_esd(alpha_esd)
                beta_esd = convert_esd(beta_esd)
                gamma_esd = convert_esd(gamma_esd)
                
                parameters.append([
                    filename, space_group_number,
                    temp_val,
                    a_val, a_esd,
                    b_val, b_esd,
                    c_val, c_esd,
                    alpha_val, alpha_esd,
                    beta_val, beta_esd,
                    gamma_val, gamma_esd,
                    vol_val,
                    Rone
                ])
    return parameters

def select_folder():
    folder_path = filedialog.askdirectory()
    if folder_path:
        entry.delete(0, tk.END)
        entry.insert(tk.END, folder_path)

def process_files():
    folder_path = entry.get()
    if folder_path:
        parameters = extract_parameters(folder_path)
        output_filename = "parameters.csv"
        output_filepath = os.path.join(folder_path, output_filename)
        with open(output_filepath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "Filename", "Space Group Number",
                "Temp",
                "a", "a (ESD)",
                "b", "b (ESD)",
                "c", "c (ESD)",
                "alpha", "alpha (ESD)",
                "beta", "beta (ESD)",
                "gamma", "gamma (ESD)",
                "Volume",
                "R1"
            ])
            writer.writerows(parameters)
        print(f"Parameters extracted and written to {output_filepath}")

root = tk.Tk()

label = tk.Label(root, text="Select a folder:")
label.pack()

entry = tk.Entry(root, width=50)
entry.pack()

button = tk.Button(root, text="Browse", command=select_folder)
button.pack()

process_button = tk.Button(root, text="Process Files", command=process_files)
process_button.pack()

root.mainloop()