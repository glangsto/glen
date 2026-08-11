#/bin/python
#Function to stack in velocity the plots and sum the total.
#HISTORY
#26Aug11 GIL give up on plotting all lines and the sum, if too may lines
#26Aug10 GIL fix labels for different plot modes
#26Aug08 GIL stack multiple model files, clean up fitting plots
#26Aug05 GIL fix velocity type, clean up argument usage.
#26Apr26 GIL return all parameters of the fit
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
        args.molecule = ""
    else:
        args.molecule = str( args.molecule)
        
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
    # set spacing between plot lines
    offset = args.offset
    nPlot = 0
    usedMaxWeight = 0.
    usedMaxFreq = 0.

    minFreq = min(rest_freqs)
    maxFreq = max(rest_freqs)
    nFreq = len(rest_freqs)
    dFreq = maxFreq - minFreq
    # set text spaceing for plotting labels of intensity vs frequency
    dFreqText = dFreq/(nFreq*10.)

    sumSigma2 = 0.
    #total weights of spectra used in sum.
    weightSum = 0.
    titleFont = 20
    annotateFont = 16
    axisFont = 18
    axisLabelFont = 20
    yMin = 1.e9    # start with huge y minimum, and find actual minimum
    
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
            
            # if this line model is a significant contribution to the model.
            if w > usedMaxWeight:
                usedMaxWeight = w
                usedMaxFreq = nu_rest
                usedMaxLabel = labels[i]

            # now update plot label offsets for data ranges
            sliceMax = I_slice.max()
            # update the spacing for top text labels
            dText = offset/5.

            if nPlot == 0:
                offset = float(args.offset)
                offset = max(sliceMax,offset)
                peakMax = max(peakMax, 7.*sliceMax)

            # now if a single maximum is way above other features
            if args.ignore and (sliceMax > peakMax):
                print("Likely a confusing line in the sum, ignoring %10.3f" %
                      (nu_rest))
                continue

              # Convert to velocity
            freq = np.array(freqs[iObs])
            vel = freq_to_velocity( freq, nu_rest)

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
                if nPlot == 0:
                    xs = f_slice.tolist()
                    ys = I_slice.tolist()
                else:
                    xs.extend(f_slice.tolist())
                    ys.extend(I_slice.tolist())
                plt.plot( f_slice, I_slice, lw=3)
                flabel = "%d: %.3f - %s" % (iObs+1, nu_rest, labels[i])
                # label line 
                plt.text(f_slice[nf-1]+dFreqText, 0., flabel, rotation=90,
                         fontsize=annotateFont)
                
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
            nPlot = nPlot + 1

            print("Obs: %d - Line Weight %6.3f for frequency %10.3f (sigma %.3f)" % \
                  (iObs+1, w, nu_rest, sigma))

            if args.plot == 'both' or args.plot == 'stack':
                # Plot stacked spectrum
                plt.plot(v_slice, I_slice + plotOffset, lw=2)
                plt.text(args.vmin + 1., plotOffset + dText,
                         f"{nu_rest:.3f}", fontsize=annotateFont)
                flabel = "%s:" % str(iObs+1)
                plt.text(args.vmin, plotOffset + dText,
                         flabel, fontsize=annotateFont)
                plt.text(args.vmax - 3., plotOffset + dText,
                         labels[i], fontsize=annotateFont)

                if args.gauss:
                    fits, final = iterative_three_gaussian_fit(v_slice, I_slice)

                # must give rest of plots some space.
                yLastMax = max(I_slice)
                # find the minium in this data set
                iMin = min(I_slice)
                # find mimimum in all data sets)
                yMin = min( yMin, iMin)
                dOffset = max( offset, yLastMax*.33)
                plotOffset = plotOffset + dOffset

            # end for all observations
            
        # end for all rest frequencies
        
    if isinstance( args.molecule, str):
        moleculeNoDollar = args.molecule.replace("$","")
    else:
        moleculeNoDollar = ""
    moleculeNoLatex = moleculeNoDollar.replace("{","")
    moleculeNo = moleculeNoLatex.replace("}","")

    # create a file name with time and molecule, but without special characters
    datetime_iso = datetime.now().replace(microsecond=0).isoformat()

    # set default values for results
    Vpeak = 0.
    Ipeak = 0.
    Isum  = 0.
    Vfit = 0.
    Vwidth = 0.
    Widthfit = 0.
    SumRms = 0.
    Vrms = 0.
    WidthRms = 0.

    # if plotting versus frequency, no summing is possible.
    if args.plot == 'freq':
        plt.xlabel("Frequency (MHz)", fontsize=axisLabelFont)
        plt.ylabel("Intensity (K)", fontsize=axisLabelFont)
        if args.title == None:
            plt.title("Frequency vs Intensity for %s" % (args.molecule), fontsize=titleFont)
        else:
            plt.title(args.title, fontsize=titleFont)
        plt.plot( [minFreq,maxFreq],[0.,0.], 'g--', lw=.5)
        # scale the range to trim x range default
        dFreq = dFreq / 25.
        plt.xlim( minFreq-dFreq,maxFreq+dFreq)
#        plt.ylim(bottom=yMin)
        # finally show result
        plt.tight_layout()
        plt.show() 
        saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)
        plt.savefig(saveFile)
        saveFile = "%s-%s.png" % (moleculeNo, datetime_iso)
        plt.savefig(saveFile)
        # exit, no summing possible
        return Vpeak, Ipeak, Isum, Vfit, Widthfit, SumRms, Vrms, WidthRms

    # now normalize weights for strongest line
    if nPlot < 1 or sumSigma2 <= 0.:
        print("No Lines in Observation Range, Exiting!")
        print("")
        exit()
    summed = summed / sumSigma2
    rms = np.std(summed)
    # final labels are based on the sum maximum and a few sigma
    yMax = max(summed)
    # prepare to calculate y-axis values for text labels
    peakOffset = max(yMax, 7.*rms, offset)
    # the text vertical spacing should be 1/25th of total Y value
    dText = (plotOffset + peakOffset)/24.

    print("Weight Sum: %.3f.  Sum of weighted sigma squares %.3f" % (weightSum, sumSigma2 ))
    
    if args.molecule == None:
        args.molecule = ""
    else:
        args.molecule = str( args.molecule)

    saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)

    # decided to give up on plotting Sum in case of too many
    # contributed lines.
    if nPlot > args.maxplot:
        # if was plotting both lines and sum
        if args.plot == 'both':
            args.plot = 'stack'
            print("Too many lines to plot the line sum also!")
            print("Re-run glen.py with argument")
            print("     --plot sum")
            print("To view the Sum of all lines")

    # if only plotting the sum, there is no offset
    if args.plot == 'sum':
        plotOffset = 0.

    dText = (plotOffset + peakOffset)/24.

    if args.plot == 'both' or args.plot == 'sum':
        # now plot average
        plt.plot(vel_grid, summed + plotOffset, lw=3)
        # put text above plotted lines
        plt.text(args.vmax - 3., plotOffset + dText,
                 "Weighted Sum", fontsize=annotateFont)

        # plot the zero line for average
        plt.plot(vel_grid, (0.00001*summed) + plotOffset, 'g--', lw=.5)

        plt.ylabel("Weighted Intensity (K)", fontsize=axisLabelFont)
        # plot the zero line for average
        
        if args.title == None:
            plt.title("Stack in Velocity %s " % (args.molecule), fontsize=titleFont)
        else:
            plt.title(args.title, fontsize=titleFont)
    
    plt.xlim(args.vmin, args.vmax)        
    plt.xlabel("Velocity (km/s)", fontsize=axisLabelFont)
    plt.ylabel("Intensity (K) + offset", fontsize=axisLabelFont)
    # label for top of plot default
    plt.ylim(bottom=-offset/2.)
    
    if args.title == None:
        plt.title("Stack in Velocity %s " % (args.molecule), fontsize=titleFont)
    else:
        plt.title(args.title, fontsize=titleFont)
    
        
    # Now fit a gaussian to the sum
    sumMax = summed.max()
    iMax = summed.argmax()
    velMax = vel_grid[iMax]

    Widthfit = 0.
    SumRms = 0.
    Vrms = 0.
    WidthRms = 0.
    vrange = args.vmax - args.vmin

    # initial guess is at max location, with 1/20th toe velocity rangea
    initial_guess = [ sumMax, velMax, ((vrange)*.05)]
    # set default search velocity
    if args.velocity != None:
        velocity = float( arg.velocity)
    else:
        velocity = 0.


    fit_A = 0.
    # now try to fit sum velocity
    try:
        popt, pcov = curve_fit(gaussian, vel_grid, summed, p0=initial_guess)
        # Extract the fitted parameters
        fit_A, fit_x0, fit_sigma = popt
        errors = np.sqrt(np.diag(pcov))
        print("Fitted parameters: Peak %10.6f +/- %10.6f" % (fit_A, errors[0]))
        print("                 : Velocity   %6.3f km/sec" % (fit_x0))
        print("                 : Line Width %6.3f km/sec" % (fit_sigma))
        # if velocity not specified and fit seems successful
        wlabel = ""
        if args.velocity == None:
            vlabel = "Vel: %.3f$\pm$%.3f" % \
                (fit_x0, errors[1])
            wlabel = "Width:%.3f$\pm$%.3f" % \
                (fit_sigma, errors[2])
            Vpeak = fit_x0
            WidthFit = fit_sigma
            velocity = Vpeak
        else:
            Vpeak = float(args.velocity)
            vlabel = "Vel: %.3f" % (float(args.velocity))
            velocity = Vpeak

        # sanity check on fit.  If peak is a small sigma, just use max.
        if fit_A < 2.5 * rms:
            print( "Peak fit is not significant, reporting max instead")
            iMax = np.argmax( summed)
            iMin = np.argmin( summed)
            if summed[iMax] > - summed[iMin]:
                fit_A = summed[iMax]
                Vpeak = vel_grid[iMax]
            else:
                fit_A = summed[iMin]
                Vpeak = vel_grid[iMin]
            plabel = "Peak: %.3f$\pm$%.3f" % (fit_A, rms)
            vlabel = "Vel: %.3f" % (Vpeak)
            fitOK = False
        else:        
            print("Max Model Line   : %12.6f   %s (%.3f)" %
                  (usedMaxFreq, usedMaxLabel, usedMaxWeight))
            plabel = "Peak: %.3f$\pm$%.3f" % (fit_A, errors[0])
            fitOK = True

        plt.axvline(Vpeak, color='blue', linestyle="--", linewidth=1)

        topMax = max( 7.*rms, sumMax, offset)
        dText = (plotOffset+topMax)/24.
        # trying to place sum labels in correct locations.
        if args.plot == 'sum':
            plotOffset = 0.

        # need to estimate the Y axis range to put text in correct spot
        if args.plot == 'both' or args.plot == 'sum':
            Vtext = Vpeak - (vrange*.45)
            if Vtext < args.vmin:
                Vtext = args.vmin + (vrange*.1)
            plt.text(Vtext, plotOffset+1.*dText,
                     vlabel, fontsize=annotateFont)
            plt.text(Vtext, plotOffset+2.*dText,
                     wlabel, fontsize=annotateFont)

            plt.text(Vtext, plotOffset+3.*dText,
                     plabel, fontsize=annotateFont)
            if fitOK:
                plt.plot(vel_grid, gaussian(vel_grid, *popt)+plotOffset, 'r--')
            # now do not double plot velocity
        Vpeak = None
        
    except RuntimeError as e:
        print(f"Error during fitting: {e}. No summed intensity detected.")

    if args.velocity != None:
        Vpeak = float(args.velocity)
        vlabel = "Vel: %.2f" % (Vpeak)
        plt.axvline(velocity, color='blue', linestyle="--", linewidth=1)
        plt.text(velocity-(vrange*.4), plotOffset+sumMax-dText,
                 vlabel, fontsize=annotateFont)
    else:
       Vpeak = 0.0
    # If just stacking, just need space for last plot
    if args.plot == 'stack':
        peakOffset = yLastMax
        
    if args.plot == 'both':
        plt.ylim(bottom=-offset/3.,top=plotOffset+peakOffset+dText)
    else:
        plt.ylim(bottom=-offset/3.,top=plotOffset+peakOffset)

    plt.tight_layout()
    # finally show result
    plt.show()

    saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    saveFile = "%s-%s.png" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    
    # Report S/N improvement
    if args.report_snr:
        rms_individual = np.mean(individual_rms)
        rms_sum = np.std(summed)
        theoretical_gain = np.sqrt(nPlot)
        measured_gain = rms_individual / rms_sum if rms_sum > 0 else np.inf
        rms1 = np.std(summed[0:100])
        rms2 = np.std(summed[(nSum-100):(nSum-1)])
        rms = max(rms1,rms2)

        print("\n--- S/N Improvement --- %10.5f" % (measured_gain))
        print("RMS of sum %10.5f (%10.5f, %10.5f)" % (rms, rms1, rms2))

    return Vpeak, Ipeak, Isum, Vfit, Widthfit, SumRms, Vrms, WidthRms
