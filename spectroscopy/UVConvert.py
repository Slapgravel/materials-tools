import os
import csv
from specio import specread
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# Specify the folder containing your spectra files
folder_path = r"\\snc-isiarchive03-smb.thefacebook.com\xroptics_labshare_smb_001\Perkin_Elmer\Jason_OSC\2.7OB"

min_wavelength = 380 # adjust this value based on your needs
max_wavelength = 700 # adjust this value based on your needs



# Iterate through all files in the folder with the .sp extension
for filename in os.listdir(folder_path):
    if filename.endswith(".sp"):
        # Construct the full path to the file
        file_path = os.path.join(folder_path, filename)
        
        # Read the spectra data from the file
        spectra = specread(file_path)

        csv_filename = filename.replace(".sp", ".csv")
        csv_path = os.path.join(folder_path, csv_filename)
        with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)

            # Write the metadata as a header
            for key, value in spectra.meta.items():
                value_str = str(value).replace('\x00', '')
                writer.writerow([f"# {key}: {value_str}"])

            # Write the header row (wavelength and amplitude)
            writer.writerow(["Wavelength", "Amplitude"])

            # Write each data point as a row in the CSV file using the entire data range
            for i in range(len(spectra.wavelength)):
                writer.writerow([spectra.wavelength[i], spectra.amplitudes[i]])

        # Filter the spectra data
        idx = np.where((spectra.wavelength >= min_wavelength) & 
                      (spectra.wavelength <= max_wavelength))[0]
        sliced_wavelength = spectra.wavelength[idx]
        sliced_amplitudes = spectra.amplitudes[idx]
        
        # Make sure your data is properly scaled
        # If your data is in %R (0-100), convert to decimal (0-1)
        if np.max(sliced_amplitudes) > 1:
            sliced_amplitudes_decimal = sliced_amplitudes / 100
        else:
            sliced_amplitudes_decimal = sliced_amplitudes
            
        # Calculate absorption and Kubelka-Munk transform
        absorption = 1 - sliced_amplitudes_decimal
        F_R = (1 - absorption)**2 / (2 * absorption)
        E = 1240 / sliced_wavelength
        fn_trunc = filename.split('.', 3)
        print(fn_trunc)
        
        # REPLACE THIS SECTION with the new delta-finding code:
        if len(F_R) > 20:  # Ensure there's enough data points
            # Find the section of the data with the largest delta (over a window of 10 points)
            max_delta_idx = 0
            max_delta = 0
            window_size = 10
            
            for i in range(len(F_R) - window_size):
                delta = F_R[i+window_size] - F_R[i]
                if abs(delta) > max_delta:
                    max_delta = abs(delta)
                    max_delta_idx = i
            
            # Fit a straight line to the data points around the section with the largest delta
            X = E[max_delta_idx:max_delta_idx+window_size].reshape(-1, 1)
            y = F_R[max_delta_idx:max_delta_idx+window_size]
            
            model = LinearRegression()
            model.fit(X, y)
            
            # Calculate the x-intercept of the trendline (bandgap energy)
            bandgap = -model.intercept_ / model.coef_[0]
            
            # Plot the original spectra
            fig, axs = plt.subplots(2, figsize=(8, 6))
            axs[0].plot(sliced_wavelength, sliced_amplitudes)
            axs[0].set_xlabel('Wavelength (nm)')
            axs[0].set_ylabel('Reflectance')
            axs[0].set_title(f'{fn_trunc[0]} ({fn_trunc[2]})')
            # Plot the Kubelka-Munk transformed spectra
            axs[1].plot(E, F_R)
            axs[1].set_xlabel('Energy (eV)')
            axs[1].set_ylabel('Kubelka-Munk Function')
            axs[1].set_title('Kubelka-Munk Transformed Spectra')
            axs[1].set_ylim(-1, 10)
            axs[1].axhline(0, color='black', linestyle='--')
            
            # Plot the fitted line - ensuring it covers the right range
            x_range = np.linspace(min(E), max(E), 100).reshape(-1, 1)
            axs[1].plot(x_range, model.predict(x_range), "r--")
            
            # Highlight the points used for the trendline
            axs[1].scatter(E[max_delta_idx:max_delta_idx+window_size], 
                          F_R[max_delta_idx:max_delta_idx+window_size], 
                          color='green', s=30)
            
            # Print the bandgap
            print(f"Bandgap: {bandgap:.2f} eV")
            
            # Add the bandgap value to the plot
            axs[1].text(0.95, 0.95, f"Bandgap: {bandgap:.2f} eV", 
                       transform=axs[1].transAxes, fontsize=10,
                       verticalalignment='top', horizontalalignment='right')
            
            # Layout so plots do not overlap
            fig.tight_layout()
            
            # Show the plot
            plt.show()
        else:
            # Code for when there aren't enough data points
            fig, axs = plt.subplots(2, figsize=(8, 6))
            axs[0].plot(sliced_wavelength, sliced_amplitudes)
            axs[0].set_xlabel('Wavelength (nm)')
            axs[0].set_ylabel('Reflectance')
            axs[0].set_title(f'{fn_trunc[0]} ({fn_trunc[2]})')
            
            axs[1].plot(E, F_R)
            axs[1].set_xlabel('Energy (eV)')
            axs[1].set_ylabel('Kubelka-Munk Function')
            axs[1].set_title('Kubelka-Munk Transformed Spectra')
            
            plt.show()
            print("Cannot calculate bandgap for this file. Not enough data points.")