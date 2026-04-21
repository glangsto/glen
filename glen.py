#python
#This module plots stacked versions of molecular lines
#based on Gotham and other GBT Spectra 
#HISTORY
#26Apr21 GIL merge weights and multiple spectra fitting
#26Apr11 GIL use TMC Turner-Langston data by default
#26Feb26 GIL start to use weights in Sum, add 3rd plot type which is freq vs intensity
#26Feb25 GIL Assume 1st and 2nd FITS file columns are Freq and Temp
#26Feb24 GIL ArgParse returns structre not individual items
#26Feb23 GIL check for astropy before import
#26Feb21 GIL add more command line arguments
#26Feb20 GIL use parsearg and take many input parameters
#26Feb16 GIL reduce default velocity grid
#26Feb15 GIL add summing of signficant input spectra
#26Feb15 GIL initial version

import sys
import os
import matplotlib.pyplot as plt
#frequently need an additional astropy library
try:
    from astropy.table import Table
except:
    print("'astropy' must be installed first!")
    print("try:")
    print("sudo apt-get install python3-astropy")
    print("")
    exit()
import numpy as np
from pathlib import Path
from stackSumVelocity import *
from idl2mathtext import *
from getlists import *
from parselist import *
from stackArgParse import *

def getObservationFile( testName):
    """
    getObservationFile() will try to find the (large) observation files.
    Two files are included by default
    If the file is not found, then the program will try to download it
    to the current directory.
    """
    # first see if the file is in the current directroy
    testFile=str(Path.cwd()) + "/data/" + testName
    print("Looking for %s: trying %s" % (testName, testFile))
    if os.path.exists(testFile):
        spectraFile = testFile
    else:
        testFile=str(Path.cwd()) + "/" + testName
        if os.path.exists(testFile):
            spectraFile = testFile
        else:
            testFile=str(Path.home()) + "/Downloads/" + testName
            if os.path.exists(testFile):
                spectraFile = testFile
            else:
                # else the file can be downloaded.
                import subprocess
                # else maybe glangsto downloaded it (for GBO only)
                if testName == "gotham_drv.fits":
                    webFile = "https://www.gb.nrao.edu/GbtLegacyArchive/GOTHAM/calibrated/gotham_drv.fits"
                    print("No Gotham file yet found!")
                    print("Wget-ting %s" % ( webFile))
                    print("This may take a while...")
                    print("")
                    try:
                        testFile = str(Path.cwd()) + "/" + testName
                        getCommand=["/bin/wget", "-O", testFile, webFile]
                        result=subprocess.run( getCommand)
                    except:
                        print("Failed to get a spectra measurement file")
                        exit()
                print("wget was successful, continuing")
                spectraFile=testFile
    print("%s -> %s" % (testName, spectraFile))
    return spectraFile
    # end of getObservationFile()
    
def getObservations( testNames):
    """ 
    getObservations looks for two "standard" GBT Observation files and
    if any files found, returns a list of file names
    """

    fileList = []
    nFiles = 0
    nTry = len(testNames)
    # for all file names to search file
    for iTry in range(nTry):
        testName = testNames[iTry]
        spectraFile = getObservationFile( testName)
        if spectraFile != "":
            fileList.append(spectraFile)
            nFiles = nFiles + 1
    if nFiles < 1:
        print("No Observing file found, exiting")
        exit()
    return fileList
    # end of getObservations)

# -1. Parse all the input arguments
args = stackArgParse()

spectraFiles = args.data_files.split(" ")
# if no spectra observation supplied, try to find Turner-Langston and Gotham obs, the default
spectraFiles = getObservations( spectraFiles)

nObs = len( spectraFiles)
if nObs < 1:
    print("!!! Can not find any observation files matching: %s" % (spectraFiles))
else:
    print("Attempting line stacking for %d observation files" % (nObs))
    
# 0. First read a line list
# The line list is formated for use in idl.  This list will be converted
# into a matplotlib friendly array

if args.intensity == None:
    print("A model intensity file is needed")
    args.intensity = "pro/hc5n.pro"
    args.molecule = "$HC_5N$"
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

names = []
nGood = 0   # count number of found observations
#devide the files by spaces
    
for iObs in range(nObs):
    # 1. Read the FITS file into an Astropy Table object
    # The Table.read() method is a high-level interface for reading FITS tables.
    # prepare to read several observing files

    try:
        table = Table.read(spectraFiles[iObs], format='fits')
    except FileNotFoundError:
        print("Error: The file %s was not found." % (spectraFiles[iObs]))
        exit()
    except Exception as e:
        print(f"Error reading the FITS file: {e}")
        exit()

    # for all spectral data files, abbreviate known Surveys
    if "gotham" in spectraFiles[iObs].lower():
        names.append("GOTHAM")
    else:
        if "tmc-tlq" in spectraFiles[iObs].lower():
            names.append("TMC-TLQ")
        else: # deduce survey name from file name
            # else unknow survey, parse file name
            aFile = spectraFiles[iObs]
            aparts = aFile.split("/")
            nparts = len( aparts)
            filepart = aparts[nparts-1]
            aparts = filepart.split(".")
            surveyName = aparts[0]
            names.append(surveyName)

    # sanity test, try to find object name
    try:
        print("Object: %s" % (table.meta['OBJECT']))
    except:
        print("Object Name not found")

    print("Survey %2d: File: %s" % (iObs, spectraFiles[iObs]))

    try:
        findex = table.colnames.index("frequency")
        fname  = "frequency"
    except:
        try:
            findex = table.colnames.index("FREQUENCY")
            fname  = "FREQUENCY"
        except:
            #Frequency column not found assume 1nd column
            findex = 0
            fname = table.colnames[findex]

    # try to find T main beam column. 
    try:
        tindex = table.colnames.index("Tmb")
        tname = "Tmb"
    except:
        try:
            tindex = table.colnames.index("TMB")
            tname = "TMB"
        except:
            #Temp column not found assume 2nd column
            tindex = 1
    tname = table.colnames[tindex]
    # now use indices to get column units
    xaxis=table.columns[fname].unit
    yaxis=table.columns[tname].unit
    
    # 2. Access the data columns
    # You can access columns by their names.
    # Replace 'column_name_x' and 'column_name_y' with the actual column names in your file.
    try:
        x_data = table[fname]
        y_data = table[tname]

        # Convert to plain NumPy arrays if needed for specific plotting functions
        x_array = np.array(x_data)
        y_array = np.array(y_data)
        nGood = nGood + 1
    except KeyError as e:
        print("For observation %2d: Can not parse data in file %s" % (iObs, spectraFiles[iObs]))
        print(f"Error: Column name {e} not found in the FITS table.")
        print(f"Available columns are: {table.colnames}")
        continue
        
# now loop and stack velocities for multiople line lists
    nCol = len( table.colnames)    # get count of columns
    freqName = table.colnames[findex]   # this assumes first column is Frequency
    tempName = table.colnames[tindex]   # second column is intensity
        
    print("Plotting %s vs %s " % (tempName, freqName))
    if nCol > 2:                   # if an RMS/sigma column
        rmsName = table.colnames[2]

        rms_array = table[rmsName]
        nRms = len( rms_array)
        nRms2 = int(nRms/2)
        rms = [rms_array[nRms2]]
    else:
        # if no rms column compute RMS
        nY = len( yaxis)
        nbY = int(nY/10)
        neY = nY - nby
        rms = np.stddev( yaxis[nbY:neY])

    # if first observing file, create the data lists
    if iObs == 0:
        xs = [x_array]
        ys = [y_array]
        rmss = [rms]
        # only will need one representative RMS for each sum.
    else:  # else append to lists
        xs = [xs, x_array]
        ys = [ys, y_array]
        rmss = [rmss, rms]
        
# end of reading all observing files

if nGood < 1:
    print("!!! Unable to read any spectra, exiting !!!")
    exit()
    
# convert labels to matplot lib nice formats
for i, alabel in enumerate( labels):
    labels[i] = idl2mathtext( alabel)

if args.title == None:
    if args.plot == 'freq':
        args.title = "%s Line Strength Comparison for %s" % (args.molecule, str(names))
    else:
        args.title = "Stack %s in Velocity %s Observations" % (args.molecule, str(names))

# finally do all computations
velocity, intensity = stackSumVelocity( xs, ys, rmss, nObs,
                                        freqs, weights, labels, args)

