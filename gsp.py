import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

# 1. Read the FITS file into an Astropy Table object
# The Table.read() method is a high-level interface for reading FITS tables.
# Replace 'your_fits_file.fits' with the path to your file.
try:
    gothamFile='/users/glangsto/Downloads/gotham_drv.fits'
    table = Table.read(gothamFile, format='fits')
    print(f"Successfully read table with columns: {table.colnames}")
except FileNotFoundError:
    print("Error: The file %s was not found." % (gothamFile))
    exit()
except Exception as e:
    print(f"An error occurred while reading the FITS file: {e}")
    exit()

print(table.meta)
print(table.columns)
print(table.columns['frequency'].unit)
xaxis=table.columns['frequency'].unit
yaxis=table.columns['Tmb'].unit
    
# 2. Access the data columns
# You can access columns by their names.
# Replace 'column_name_x' and 'column_name_y' with the actual column names in your file.
try:
    x_data = table['frequency']
    y_data = table['Tmb']

    # Convert to plain NumPy arrays if needed for specific plotting functions
    x_array = np.array(x_data)
    y_array = np.array(y_data)

except KeyError as e:
    print(f"Error: Column name {e} not found in the FITS table.")
    print(f"Available columns are: {table.colnames}")
    exit()

# 3. Plot the data using Matplotlib
plt.figure(figsize=(8, 6))
#plt.scatter(x_array, y_array, marker='.')
plt.plot( x_array, y_array)
plt.xlabel(xaxis)
plt.ylabel(yaxis)
plt.title('Scatter Plot of Gotham GBT survey of TMC-1')
plt.grid(True)
plt.show()
