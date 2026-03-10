import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Read the CSV file into a Pandas dataframe
df = pd.read_csv('data.csv', header=0)              # Set header line number

# Get the column names
column_names = df.columns
# df['diff'] = df[column_names[2]] - df[column_names[1]]*0.1*df[column_names[3]]+0.5                # Random manipulation of data
# df['dydx'] = np.gradient(df[column_names[2]])           # Differentiate data set
fig, ax = plt.subplots()

# Plot the data
plt.plot(df[column_names[0]], df[column_names[1]], label=column_names[1])
plt.plot(df[column_names[0]], df[column_names[2]], label=column_names[2])
# plt.plot(df[column_names[0]], df['dydx'], label='dydx')
plt.xlabel(column_names[0])
plt.ylabel(column_names[1])
ax.axhline(y=0, color='black', linewidth=0.5, linestyle='--')               # horizontal line @ y=0
ax.set_facecolor((0.17, 0.17, 0.17))                # plot backgound colour
fig.patch.set_facecolor((0.15, 0.15, 0.15))             # fig background color
ax.set_title("My Plot", fontdict={'family': 'sans-serif', 'color': 'white', 'size': 14})              # Set the title with a specific font and color
ax.set_xlabel(r"Degrees ($^\circ$)", fontdict={'family': 'sans-serif', 'color': 'white', 'size': 14})             # Set the x and y labels with a different font and color
ax.set_ylabel("Y-axis", fontdict={'family': 'sans-serif', 'color': 'white', 'size': 14})
ax.tick_params(axis='both', labelsize=12, labelcolor='white', labelfontfamily='sans-serif')             # Set tick label style
legend = ax.legend(shadow=True, fontsize='small')             # Set legend
plt.grid(axis='both', linestyle='--', linewidth=1, color=(0.25, 0.25, 0.25))             # gridlines
plt.grid(True)
# plt.fill_between(df[column_names[0]], df[column_names[1]], df[column_names[2]],alpha=0.2,color='black')
# plt.fill_between(df[column_names[0]], df['dydx'], df[column_names[3]],alpha=0.2,color='blue')
plt.show()