#/bin/python
#Model molecular abundance based on GBT observations of molecular line
#This script includes all the assumptions needed to convert to mass
#This script is designed for linear rotator calculations of gass phase
#Molecules
# intensity

#HISTORY
#26Apr26 GIL Document test case
#26Apr23 GIL Initial version

import numpy as np
from scipy.optimize import curve_fit

# masses for for comparision
Mmilkyway = 2.85e45   #Milky Way mass grams
Msun   = 1.988475e33  #Suns mass in grams 
Mearth = 5.97e27      # Earths mass grams
Mmoon  = 7.3477       # Earth Moon Mass grams.
Mocean = 1.4e24       # Earth Ocean Mass grams.
#PlanckHBoltzmannK = 4.799e-17 #  h/k MHz/K
BoltzmannK = 1.381e-16          #  ergs/K
PlanckHBoltzmannK = 20836.6     #  MHz/K
PlanckBoltzmann = 4.798e-11      #  cgs, K seconds
Nfactor = 1.737e15    # GHz^2 Debye^2

def beamsize( frequencyMHz):
    """
    Model the GBT beamsize for the observing frequency in MHz
    returns the beam size in arcseconds
    """
    if frequencyMHz <= 0.:
        print(" Invalid input frequency: %.3f" % (frequencyMHz))
        return 3600.  # guess a degree
    return 720. * 1000. / frequencyMHz


def omegaB( frequencyMHz):
    """ 
    Beam area steradians from GBT beam size arc seconds
    """
    beamAs = beamsize( frequencyMHz)
    beamRad = beamAs * np.pi / (180.*3600.)
    
    return 1.133 * beamRad * beamRad

def distanceCm( distanceLy):
    """
    Convert distance in Light years to distance in cm
    """
    return distanceLy * 9.461e+17

def beamArea( frequencyMHz, distanceLy):
    """
    Compute GBT beam area from frequency and distance
    """
    return omegaB( frequencyMHz) * distanceCm( distanceLy)**2

def integratedIntensity( peakK, widthKmSec):
    """
    Compute the integrated intensity from gaussian fit to line
    """
    return peakK * widthKmSec * np.sqrt(np.pi/(4. * np.log(2.)))

def gUp( Jup):
    """
    Molecule partition function calculation
    """
    return (2.*Jup) + 1.

def S(Jup):
    """
    Molecule S function calculation
    """
    return Jup / ((2.*Jup) + 1)

def EuK( Jup, BfrequencyMHz):
    """
    Upper Energy Level, based on quantum number and interval between lines
    returns Eu/k in Kelvins
    """
    return PlanckBoltzmann * BfrequencyMHz * 1.e6 * Jup * (Jup + 1)

def Q( Tphysical, BfrequencyMHz):
    """ 
    Partition function for a linear rotator at Low temperature
    """
    return Tphysical/(PlanckBoltzmann * BfrequencyMHz * 1.e6)

def Ncm2( freqLineMHz, Tint, BfrequencyMHz, Tphysical, Jup, DipoleDebye, verbose=False):
    """
    Calculate the total number of molecules per cm**2
    from the measured parameters
    This approximation is appropriate for optically thin, low energy
    emission in thermal equilibrium.
    where
    freqLineMHz   is the measured line frequency in MHz
    Tint          is the integrated line intensity in K km/sec
    BfrequencyMHz is the linear component of spacing between molecular lines
    Tphysical     is the physical temperature of the emitting region
    Jup           is the upper quantum number of the transition from Jup to Jup-1
    DipoleDebye   is the molecular dipole moment
    """
    DebyeCgs = DipoleDebye * 1.e-18
#   part = (3./PlanckHBoltzmannK) / (8. * np.pi * np.pi * np.pi * freqLineMHz * S( Jup) * DebyeCgs**2)
    part = Nfactor/(BfrequencyMHz*DipoleDebye/1.E3)**2
#   From institute of Astronomy:  
    expPart = np.exp( 0.048*(BfrequencyMHz/1e3)*(Jup)*(Jup+1.)/Tphysical)
    Eu = EuK( Jup, BfrequencyMHz)
    if verbose:
        print( "Jup: %.2f S( Jup): %.2f " % (Jup, S(Jup)))
        print( "Eu/k: %.3e, Eu/kT: %.3e" % (Eu, Eu/Tphysical))
#    Ntot = part * Q(Tphysical, BfrequencyMHz) * np.exp( Eu/Tphysical) * Tint
    Tplus = Tphysical + (0.016*BfrequencyMHz/1.e3)
    if verbose:
        print("part: %.3e  expPart: %.3e  Tplus: %.3e " % (part, expPart, Tplus))
    
    Ntot = part * (expPart / Jup**2) * Tplus * Tint

    if verbose:
        print( "Q: %.3f, Exp: %.3f, Tint: %f" % ( Q(Tphysical, BfrequencyMHz),
                                                  np.exp( EuK( Jup, BfrequencyMHz)/Tphysical),
                                                 Tint))
    # end of Ncm2
    return Ntot
    
if __name__ == "__main__":
    """ 
    Default test is use of 
    Example for HC5N Gbt Observations
    """
    object = "TMC-CP"
    dashes = "-----------------------------------------------------------------------------"
    print(dashes)
    print("Running abundance model with GBT observations of %s" % (object))
    freqLineMHz = 39939.591  #MHz
    Tpeak = 3.510  # K
    deltaV = 0.203 # km/sec
    Tint = integratedIntensity( Tpeak, deltaV)
    BfrequencyMHz = 1331.3 # MHz
    Tphysical = 5. # Kelvins
    Jup  = 15
    DipoleDebye = 4.33
    Matom = 1.22e-22 # gram

    molecule = "HC5N"
    print("Langston and Turner observations of molecule %s, line Jup: %d -> %d" % (molecule, Jup, Jup-1))
    print(dashes)
    print("Laboratory Molecule Measurements:")
    print("Molecule mass %.3e gm" % (Matom))
    print("Molecule Dipole Moment: %.3f Debyle, Rotation B Frequency: %.3f MHz" % (DipoleDebye, BfrequencyMHz))
    print(dashes)
    print("GBT Measurements:")
    print("Frequency %.3f MHz" % (freqLineMHz))
    print("Gaussian Fit Peak Intensity: %.3f K and Velocity Width %.3f km/sec" % (Tpeak, deltaV))
    print("Integrated Intensity       : %.3f K km/sec.  Model Physical Temp %.3f K" % (Tint, Tphysical))
    print(dashes)
    print("Assumptions:")
    print("Linear rotator molecule model.")
    print("Emission region is optically thin.")
    print("Entire source is within the telescope beam.")
    print("Telescope beam size model based on GBT Measurements vs Frequency")

    # distance reference:    Galli, P. A. B., Loinard, L., Bouy, H., et al. 2019, A&A, 630, A137
    distancePc = 140.2 # +/- 1.3 parsecs
    print( "Distance to %14s distance: %.3f pc" % (object, distancePc))
    print( "Distance reference: Galli, P. A. B., et al. 2019, A&A, 630, A137")

    Ntot = Ncm2( freqLineMHz, Tint, BfrequencyMHz, Tphysical, Jup, DipoleDebye)

    distanceLy = distancePc * 3.26156
    Area = beamArea( freqLineMHz, distanceLy)
    
    Nmolecule = Ntot * Area

    print(dashes)
    print( "Results:")
    print( "Number of Molecules per square cm: %.3e" % (Ntot))
    print( "GBT Beam Area at Object distance : %.3e cm^2;  Nmolecule = %.2e" % (Area, Nmolecule))
    Mtotal = Matom*Nmolecule
    print( "Mass total: %.3e (gm).  Relative to Earth's Ocean = %.3e (Oceans)" % (Mtotal, Mtotal/Mocean))
    print(dashes)
    print("Glen Langston, 2026 April 26 ")
    print(dashes)
