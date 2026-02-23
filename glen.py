#python
#This module plots stacked versions of molecular lines
#based on Gotham and other GBT Spectra 
#HISTORY
#26Feb23 GIL check for astropy before import
#26Feb21 GIL add more command line arguments
#26Feb20 GIL use parsearg and take many input parameters
#26Feb16 GIL reduce default velocity grid
#26Feb15 GIL add summing of signficant input spectra
#26Feb15 GIL initial version

import sys
import os
import matplotlib.pyplot as plt
try:
    from astropy.table import Table
except:
    print("'astropy' must be installed first!")
    print("try:")
    print("sudo apt-get install python3-astropy")
    print("")
    exit()
import numpy as np
from stackSumVelocity import *
from idl2mathtext import *
from getlists import *
from parselist import *
from stackArgParse import *

def getObservationFile():
    """
    getObservationFile() will try to find the (large) observation file.
    If the file is not found, then the program will try to download it
    to the current directory.
    """
    # first see if the file is in the current directroy
    testFile='./gotham_drv.fits'
    if os.path.exists(testFile):
        spectraFile = testFile
    else:
        from pathlib import Path

        testFile = str(Path.home()) + "/Downloads/gotham_drv.fits"
        print("TestFile: %s" % (testFile))
        # save the download file, in case the  next step does not work
        downloadFile = testFile
# else see if the file file was downloaded
        if os.path.exists(testFile):
            spectraFile = testFile
        else:
            # else maybe glangsto downloaded it (for GBO only)
            testFile='/users/glangsto/Downloads/gotham_drv.fits'
            if os.path.exists(testFile):
                spectraFile = testFile
            else:
                # else the file can be downloaded.
                import subprocess
                webFile = "https://www.gb.nrao.edu/GbtLegacyArchive/GOTHAM/calibrated/gotham_drv.fits"
                print("No Gotham file yet found!")
                print("Wget-ting %s" % ( webFile))
                print("This may take a while...")
                print("")
                try:
#                    wgetArgs = "-O %s %s" % (downloadFile, webFile)
                    getCommand=["/bin/wget", "-O", downloadFile, webFile]
                    result=subprocess.run( getCommand)
                except:
                    print("Failed to get a spectra measurement file")
                    exit()
                print("wget was successful, continuing")
                spectraFile=downloadFile
    return spectraFile
    # end of getObservationFile()
    
# -1. Parse all the input arguments
args = stackArgParse()

spectraFile = ""
# if no spectra observation supplied, try to find Gotham obs, the default
if args.data_files == None:
    spectraFile = getObservationFile()
else:
    spectraFile = args.data_files

print("Stacking spectral lines from Model file: %s" % (spectraFile))

# 0. First read a line list
# The line list is formated for use in idl.  This list will be converted
# into a matplotlib friendly array

if args.intensity == None:
    print("A model intensity file is needed")
    args.intensity = "pro/hc5n.pro"
    print("Using model default: %s" % (args.intensity))

print("Parsing Model Intensity File (IDL): %s" % (args.intensity))
with open(args.intensity) as f:
    text = f.read()

freqs, labels, weights = extract_idl_arrays(text)

print("%3d freqs   1st: %.3f" % (len(freqs), freqs[0]))
print("%3d labels  1st: %s" %  (len(labels), labels[0]))
print("%3d weights 1st: %.3f" % (len(weights), weights[0]))
f.close()

molecule= args.molecule
# if no user supplied molecue name
if molecule == None:
    try:
        # try to find specific molecule tag from model.
        declared, funcs, referenced = parse_idl_file(args.intensity)
        iPlot = declared.index('plotTitle')
        moleculeIdl = referenced[iPlot]
        molecule = idl2mathtext(moleculeIdl)
    except:
        iPlot = 0
    
# 1. Read the FITS file into an Astropy Table object
# The Table.read() method is a high-level interface for reading FITS tables.
# Replace 'your_fits_file.fits' with the path to your file.
try:
    table = Table.read(spectraFile, format='fits')
    print(f"Successfully read spectral table with columns: {table.colnames}")
except FileNotFoundError:
    print("Error: The file %s was not found." % (spectraFile))
    exit()
except Exception as e:
    print(f"An error occurred while reading the FITS file: {e}")
    exit()

try:
    print("Object: %s" % (table.meta['OBJECT']))
except:
    print("Object Name not found")
#print(table.columns['frequency'].unit)
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

# now loop and stack velocities for multiople line lists

# convert labels to matplot lib nice formats
for i, alabel in enumerate( labels):
    labels[i] = idl2mathtext( alabel)

stackSumVelocity( x_array, y_array, freqs, labels,
                  molecule, 
                  vmin=args.vmin, vmax=args.vmax, dv=args.dv,
                  offset=args.offset,
                  weights=weights, normalize=False,
                  velocity=args.velocity,
                  plot=args.plot,
                  baseline=args.baseline,
                  survey=args.survey)

