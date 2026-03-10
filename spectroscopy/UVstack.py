import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore, QtGui
import numpy as np
import csv
import tkinter as tk
from tkinter import filedialog

# Function to select files
def select_files():
    root = tk.Tk()
    root.withdraw()  # Hide the Tkinter window
    files = filedialog.askopenfilenames(title="Select CSV files", filetypes=[("CSV files", "*.csv")])
    return files

# Select files
files = select_files()

# Initialize empty lists to store data
wavelengths = []
amplitudes = []
names = []

# Read selected files
for file in files:
    name = file.split('/')[-1].split('.')[0]
    names.append(name)
    
    # Initialize empty lists for this file's data
    wavelengths_file = []
    amplitudes_file = []
    
    with open(file, 'r') as f:
        reader = csv.reader(f)
        
        # Skip header lines
        for row in reader:
            if row[0] == 'Wavelength':
                break
        
        # Read data lines
        for row in reader:
            wavelengths_file.append(float(row[0]))
            amplitudes_file.append(float(row[1]))
    
    wavelengths.append(wavelengths_file)
    amplitudes.append(amplitudes_file)

# Create the plot
app = QtWidgets.QApplication([])
win = pg.GraphicsLayoutWidget(show=True, title="Visible Spectrum")
plot = win.addPlot()

# Add visible light spectrum as a background
visible_spectrum = np.linspace(380, 750, 1000)
colors = [(255, 0, 255), (0, 0, 255), (0, 255, 0), (255, 255, 0), (255, 0, 0)]
for i in range(len(visible_spectrum) - 1):
    x = [visible_spectrum[i], visible_spectrum[i+1]]
    y = [0, 0]
    bg = pg.PlotDataItem(x=x, y=y, pen=None, brush=pg.mkBrush(colors[int((i / len(visible_spectrum)) * len(colors))]))
    plot.addItem(bg)

# Plot all files
for i in range(len(names)):
    curve = pg.PlotDataItem(x=wavelengths[i], y=amplitudes[i], name=names[i])
    plot.addItem(curve)

plot.setLabel('bottom', 'Wavelength (nm)')
plot.setLabel('left', 'Amplitude')
plot.showGrid(x=True, y=True)

# Start the application
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        app.exec()