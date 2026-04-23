#/bin/python
#Function to stack in velocity the plots and sum the total.
#HISTORY
#26Apr22 GIL bigger font for dual plots
#26Apr21 GIL clean up intensity vs frequency plot
#26Apr11 GIL fit up to 3 gaussians to input individual spectra
#26Feb26 GIL revised sums to use RMSs for optimum measurement
#26Feb25 GIL take lists of arrays of observations and look for matches 
#26Feb23 GIL just pass args structure, not individual items
#26Feb21 GIL add median baseline 
#26Feb20 GIL add fit labeling
#26Feb16 GIL reduce printouts, clean up fitting
#26Feb15 GIL initial version

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from datetime import datetime
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from threefit import *

C = 299792.458  # km/s

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def freq_to_velocity(nu_obs, nu_rest):
    """Convert observed frequencies to velocities (radio definition)."""
    return C * (nu_rest - nu_obs) / nu_rest

def checkArgs( args):
    """ Check the input arguments and update default values
    """
    # first check arguments
    if args.vmin == None:
        args.vmin = 0.0
    else:
        args.vmin = float( args.vmin)
    if args.vmax == None:
        args.vmax=12.0
    else:
        args.vmax = float( args.vmax)
        
    if args.dv == None:
        args.dv = 0.05
    else:
        args.dv = float (args.dv)
        
    if args.offset == None:
        args.offset = 0.5,
    else:
        args.offset = float (args.offset)
    
    if args.normalize == None:
        args.normalize=False
    if args.report_snr == None:
        args.report_snr = False

    if args.plot == None:
        args.plot='both'   # options are 'both', 'stack', 'freq' and 'sum'
    else:
        args.plot = args.plot.lower()

    if args.molecule == None:
        molecule = ""
    else:
        molecule = str( args.molecule)
    if args.survey == None:
        args.survey = "GOTHAM",
    else:
        args.survey = str( args.survey)
        
    if args.baseline == None:
        args.baseline = False
    else:
        args.baseline = bool( args.baseline)

    # end of checkArgs()
    return args

# the input arrays, freq_uency and intensity are large
# Gotham Arrays, with irregular frequency samples.
def stackSumVelocity( freqs, intensitys, rmss, nObs, rest_freqs, weights, labels, args):
    """
    Convert spectra to velocity space for each rest frequency,
    stack them visually, and compute a summed spectrum on a
    uniform velocity grid.

    Parameters
    ----------
    freqs : lists of frequency arrays
        Observed frequency array.
    intensitys : lists of intensity arrays
        Observed intensity array.
    rmss : list
        A single RMS for each pair of frequencies and intensities
    nObs : integer
        Number of different observation files plotting
    rest_freqs : list
        Molecular line rest frequencies.
    weights : list
        Molecular line weights for sums
    labels : list
        List of names of spectral lines

    Returns
    -------
    vel_grid : array
        Uniform velocity grid.
    summed : array
        Summed intensity on the uniform grid.
    """

    n_lines = len(rest_freqs)

    # set default arguments appropriate for TMC 1 if no user input
    args = checkArgs( args)
    
    # prepare for autoscaling
    peakMax = 0.1,

    # Default weights = 1
    if weights == None:
        weights = np.ones(n_lines)
    weights = np.array(weights)

    # normalize weights to 1 for the max weight
    maxWeight = max( weights)
    weights = weights/maxWeight
    
    # if no labels for lines, use integers
    if labels == None:
        labels = np.arange(n_lines)
        labels = str(labels)
        
    vel_grid = np.arange(args.vmin, args.vmax + args.dv, args.dv)
    summed = np.zeros_like(vel_grid)

    # Plot the summed spectrum.   if only sum, then make wider than tall
    if args.plot == 'sum':
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot()
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.xaxis.set_minor_locator(MultipleLocator(1.))
    elif args.plot == 'freq':
        fig = plt.figure(figsize=(8,10))
        ax = fig.add_subplot()
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.xaxis.set_minor_locator(MultipleLocator(1000.))
    else:
        fig = plt.figure(figsize=(8,10))
        ax = fig.add_subplot()
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.xaxis.set_minor_locator(MultipleLocator(1.))

    # init values to clean up plotting
    lastMax = 0.
    plotOffset = 0.
    offset = args.offset
    nplot = 0
    usedMaxWeight = 0.
    usedMaxFreq = 0.

    minFreq = min(rest_freqs)
    maxFreq = max(rest_freqs)
    dFreq = maxFreq - minFreq
    # set text spaceing for plotting labels of intensity vs frequency
    dFreqText = dFreq/100.

    sumSigma2 = 0.
    #total weights of spectra used in sum.
    weightSum = 0.
    titleFont = 20
    annotateFont = 12
    axisFont = 16
    axisLabelFont = 18
    
    # now for all rest frequencies, find samples
    for i, (nu_rest, w) in enumerate(zip(rest_freqs, weights)):

        # now loop through all observations
        for iObs in range( nObs):
            # Convert to velocity
            freq = np.array(freqs[iObs])
            vel = freq_to_velocity( freq, nu_rest)

            # Select velocity window
            mask = (vel >= args.vmin) & (vel <= args.vmax)
            v_slice = vel[mask]
            intensity = np.array(intensitys[iObs])
            I_slice = intensity[mask]

            #don't use measurement if no points.
            if len(I_slice) < 100:
#                print("Not enough samples for frequency %10.3f" % (nu_rest))
                continue
            
            # now update plot label offsets for data ranges
            sliceMax = I_slice.max()
            if nplot == 0:
                offset = float(args.offset)
                offset = max(sliceMax,offset)
                peakMax = max(peakMax, 20.*sliceMax)

            # now if a single maximum is way above other features
            if args.ignore and (sliceMax > peakMax):
                print("Likely a confusing line in the sum, ignoring %10.3f" %
                      (nu_rest))
                continue

            if w > usedMaxWeight:
                usedMaxWeight = w
                usedMaxFreq = nu_rest
                usedMaxLabel = labels[i]

            # Normalize if requested
            if args.normalize:
                peak = np.max(np.abs(I_slice))
                if peak > 0:
                    I_slice = I_slice / peak
                    sliceMax = sliceMax / peak

            # if removing a constant baseline, 
            if args.baseline:
                iMedian = np.median(I_slice)
                I_slice = I_slice - iMedian
                sliceMax = sliceMax - iMedan
            
            # now if plotting intensity vs frequency
            if args.plot == 'freq':
                f_slice = freq[mask]
                nf = len(f_slice)
                # these lists are only used to return total intensity vs frequency spectrum
                if nplot == 0:
                    xs = f_slice.tolist()
                    ys = I_slice.tolist()
                else:
                    xs.extend(f_slice.tolist())
                    ys.extend(I_slice.tolist())
                plt.plot( f_slice, I_slice, lw=3)
                flabel = "%d: %.3f - %s" % (iObs+1, nu_rest, labels[i])
                # label line 
                plt.text(f_slice[nf-1]+dFreqText, 0., flabel, rotation=90, fontsize=annotateFont)
                
                if args.gauss:
                    fits, final = iterative_three_gaussian_fit(f_slice, I_slice)

                
            # Interpolate onto uniform velocity grid.
            # interp requires increaseing x axis.   Frequency -> velocity flips
            interp_I = np.interp(vel_grid, np.flip(v_slice), np.flip(I_slice))

            # Weighted sum prep
            # sum is weighted by inverse sigma 2 and weights
            sigma = float(rmss[iObs]/w)
            sigma2weight = sigma*sigma
            # now sum reciprical sigma2s
            sumSigma2 = sumSigma2 + (1./sigma2weight)
            # now scale a new spectrum
            asum = interp_I * (w/sigma2weight)
            summed += asum
            # separately sum the weights used
            weightSum = weightSum + w
            # count plots to deal first update of offset
            nplot = nplot + 1

            print("Obs: %d - Line Weight %6.3f for frequency %10.3f (sigma %.3f)" % \
                  (iObs+1, w, nu_rest, sigma))

            if args.plot == 'both' or args.plot == 'stack':
                textOffset = plotOffset + (offset/4.)
                # Plot stacked spectrum
                plt.plot(v_slice, I_slice + plotOffset, lw=2)
                plt.text(args.vmin + 1., textOffset,
                         f"{nu_rest:.3f}", fontsize=annotateFont)
                plt.text(args.vmin, textOffset,
                         f"{(iObs+1):d}", fontsize=annotateFont)
                plt.text(args.vmax - 3., textOffset,
                         labels[i], fontsize=annotateFont)

                if args.gauss:
                    fits, final = iterative_three_gaussian_fit(v_slice, I_slice)

                # must give rest of plots some space.
                dOffset = max( offset, (max(I_slice)*.33))
                plotOffset = plotOffset + dOffset

            # end for all observations
            
        # end for all rest frequencies
        
    # now normalize weights for strongest line
    if nplot < 1 or sumSigma2 <= 0.:
        print("No Lines in Observation Range, Exiting!")
        print("")
        exit()
    summed = summed / sumSigma2
    if args.plot != 'freq':
        print("Weight Sum: %.3f.  Sum of weighted sigma squares %.3f" % (weightSum, sumSigma2 ))
    
    # create a file name with time and molecule, but without special characters
    datetime_iso = datetime.now().replace(microsecond=0).isoformat()
    if isinstance( args.molecule, str):
        moleculeNoDollar = args.molecule.replace("$","")
    else:
        moleculeNoDollar = ""
    moleculeNoLatex = moleculeNoDollar.replace("{","")
    moleculeNo = moleculeNoLatex.replace("}","")
    saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)

    if args.plot == 'both' or args.plot == 'sum':
        # now plot average
        plt.plot(vel_grid, summed + plotOffset, lw=3)
        # put text above plotted lines
        textOffset = plotOffset + (offset/3.)
        plt.text(args.vmax - 3., textOffset,
                 "Weighted Sum", fontsize=annotateFont)

        # plot the zero line for average
        plt.plot(vel_grid, (0.00001*summed) + plotOffset, 'g--', lw=.5)
        
    if args.plot == 'freq':
#        plt.plot(np.array(xs), np.array(ys), lw=3)
        plt.xlabel("Frequency (MHz)", fontsize=axisLabelFont)
        plt.ylabel("Intensity (K)", fontsize=axisLabelFont)
        if args.title == None:
            plt.title("Frequency vs Intensity for %s" % (args.molecule), fontsize=titleFont)
        else:
            plt.title(args.title, fontsize=titleFont)
        # finally show result
        plt.ylim(bottom=-offset/2.)
        plt.tight_layout()
        plt.show()
        saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)
        plt.savefig(saveFile)
        saveFile = "%s-%s.png" % (moleculeNo, datetime_iso)
        plt.savefig(saveFile)
        return xs, ys
    else:
        plt.xlim(args.vmin, args.vmax)        
        plt.xlabel("Velocity (km/s)", fontsize=axisLabelFont)
        if args.plot == 'sum':
            plt.ylabel("Weighted Intensity (K)", fontsize=axisLabelFont)
        else:
            plt.ylabel("Intensity (K) + offset", fontsize=axisLabelFont)
        # label for top of plot default
        if args.title == None:
            plt.title("Stack in Velocity %s " % (args.molecule), fontsize=titleFont)
        else:
            plt.title(args.title, fontsize=titleFont)
    
    # Now fit a gaussian to the sum
    sumMax = summed.max()
    iMax = summed.argmax()
    velMax = vel_grid[iMax]
    rms = np.std(summed)

    vrange = args.vmax - args.vmin
    # initial guess is at max location, with 1/20th toe velocity rangea
    initial_guess = [ sumMax, velMax, ((vrange)*.05)]
    try:
        popt, pcov = curve_fit(gaussian, vel_grid, summed, p0=initial_guess)
        # Extract the fitted parameters
        fit_A, fit_x0, fit_sigma = popt
        errors = np.sqrt(np.diag(pcov))
        print("Fitted parameters: Peak %10.6f +/- %10.6f" % (fit_A, errors[0]))
        print("                 : Velocity   %6.3f km/sec" % (fit_x0))
        print("                 : Line Width %6.3f km/sec" % (fit_sigma))
        # if velocity not specified and fit seems successful
        if args.velocity == None:
            vlabel = "Vel: %.3f$\pm$%.3f" % \
                (fit_x0, errors[1])
            wlabel = "Width:%.3f$\pm$%.3f" % \
                (fit_sigma, errors[2])
            velocity = fit_x0
        else:
            velocity = float(args.velocity)
            vlabel = "Vel: %.3f" % (args.velocity)
        print("Max Model Line   : %12.6f   %s (%.3f)" %
              (usedMaxFreq, usedMaxLabel, usedMaxWeight))
        plt.axvline(velocity, color='blue', linestyle="--", linewidth=1)
        plabel = "Peak: %.3f$\pm$%.3f" % (fit_A, errors[0])

        if args.plot == 'both' or args.plot == 'sum':
            plt.text(velocity-(vrange*.45), plotOffset+sumMax,
                     plabel, fontsize=annotateFont)
            plt.text(velocity-(vrange*.45), plotOffset+sumMax-offset,
                     vlabel, fontsize=annotateFont)
            plt.text(velocity-(vrange*.45), plotOffset+sumMax-(2.*offset),
                     wlabel, fontsize=annotateFont)
            plt.plot(vel_grid, gaussian(vel_grid, *popt)+plotOffset, 'r--')
            # now do not double plot velocity
        velocity = None
        
    except RuntimeError as e:
        print(f"Error during fitting: {e}. No summed intensity detected.")

    if args.velocity != None:
        velocity = float(velocity)
        vlabel = "Vel: %.2f" % (velocity)
        plt.axvline(velocity, color='blue', linestyle="--", linewidth=1)
        plt.text(velocity-(vrange*.4), plotOffset+sumMax-(2.*offset),
                 vlabel, fontsize=annotateFont)

    # finally show result
    plt.ylim(bottom=-offset/2.,top=plotOffset+sumMax+offset)
    plt.tight_layout()
    plt.show()

    saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    saveFile = "%s-%s.png" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    
    # Report S/N improvement
    if args.report_snr:
        rms_individual = np.mean(individual_rms)
        rms_sum = np.std(summed)
        theoretical_gain = np.sqrt(nplot)
        measured_gain = rms_individual / rms_sum if rms_sum > 0 else np.inf
        rms1 = np.std(summed[0:100])
        rms2 = np.std(summed[(nSum-100):(nSum-1)])
        rms = max(rms1,rms2)

        print("\n--- S/N Improvement --- %10.5f" % (measured_gain))
        print("RMS of sum %10.5f (%10.5f, %10.5f)" % (rms, rms1, rms2))

    return vel_grid, summed
