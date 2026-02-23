#python
#Function to stack in velocity the plots and sum the total.
#HISTORY
#26Feb21 GIL add median baseline 
#26Feb20 GIL add fit labeling
#26Feb16 GIL reduce printouts, clean up fitting
#26Feb15 GIL initial version

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

C = 299792.458  # km/s

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def freq_to_velocity(nu_obs, nu_rest):
    """Convert observed frequencies to velocities (radio definition)."""
    return C * (nu_rest - nu_obs) / nu_rest

# the input arrays, freq_uency and intensity are large
# Gotham Arrays, with irregular frequency samples.
def stackSumVelocity(
        freq, intensity, rest_freqs, labels=None,
        molecule = "",
        vmin=0.0, vmax=12.0,
        dv=0.05, offset=0.5,
        peakMax = 0.1,
        weights=None,
        normalize=False,
        report_snr=True,
        velocity=None,
        plot='both',   # options are 'both', 'stack' and 'sum'
        label=None,
        survey="GOTHAM",
        baseline=False):
    """
    Convert spectra to velocity space for each rest frequency,
    stack them visually, and compute a summed spectrum on a
    uniform velocity grid.

    Parameters
    ----------
    freq : array
        Observed frequency array.
    intensity : array
        Observed intensity array.
    rest_freqs : list
        Molecular line rest frequencies.
    labels : list
        List of names of spectral lines
    vmin, vmax : float
        Velocity window (km/s).
    dv : float
        Desired velocity resolution for the summed spectrum.
    offset : float
        Vertical spacing for stacked plot.
    weights : list or array, optional
        Weight for each line (default = equal weights).
    normalize : bool
        Normalize each spectrum before summing.
    report_snr : bool
        Compute and print S/N improvement.

    Returns
    -------
    vel_grid : array
        Uniform velocity grid.
    summed : array
        Summed intensity on the uniform grid.
    """

    n_lines = len(rest_freqs)

    # Default weights = 1
    if weights is None:
        weights = np.ones(n_lines)
    weights = np.array(weights)

    # normalize weights to 1 for the max weight
    maxWeight = max( weights)
    weights = weights/maxWeight
    
    # if no labels for lines, use integers
    if labels is None:
        labels = np.arange(n_lines)
        labels = str(labels)
        
    # Create uniform velocity grid
    vmin = float(vmin)
    vmax = float(vmax)
    dv = float(dv)
    vel_grid = np.arange(vmin, vmax + dv, dv)
    summed = np.zeros_like(vel_grid)

    # For S/N estimation
    individual_rms = []

    # Plot the summed spectrum.   if only sum, then make wider than tall
    if plot == 'sum':
        plt.figure(figsize=(8, 5))
    else:
        plt.figure(figsize=(8, 10))

    # init values to clean up plotting
    lastMax = 0.
    plotOffset = 0.
    nplot = 0
    usedMaxWeight = 0.
    usedMaxFreq = 0.
    
    #total weights of spectra used in sum.
    weightSum = 0.
    # now for all rest frequencies, find samples
    for i, (nu_rest, w) in enumerate(zip(rest_freqs, weights)):

        # Convert to velocity
        vel = freq_to_velocity(freq, nu_rest)

        # Select velocity window
        mask = (vel >= vmin) & (vel <= vmax)
        v_slice = vel[mask]
        I_slice = intensity[mask]

        #don't use measurement if no points.
        if len(I_slice) < 100:
            print("Not enough samples for frequency %10.3f" % (nu_rest))
            continue
        print("Line Weight %6.3f for frequency %10.3f" % (w, nu_rest))

        # now update off sets for data ranges
        sliceMax = I_slice.max()
        if nplot == 0:
            offset = float(offset)
            offset = max(sliceMax,offset)
            peakMax = max(peakMax, 20.*sliceMax)

        # now if a single maximum is way above other features
        if sliceMax > peakMax:
            print("Likely a confusing line in the sum, ignoring %10.3f" %
                  (nu_rest))
            continue

        if w > usedMaxWeight:
            usedMaxWeight = w
            usedMaxFreq = nu_rest
            usedMaxLabel = labels[i]
            
        # Normalize if requested
        if normalize:
            peak = np.max(np.abs(I_slice))
            if peak > 0:
                I_slice = I_slice / peak
                sliceMax = sliceMax / peak

        # if removing a constant baseline, 
        if baseline:
            iMedian = np.median(I_slice)
            I_slice = I_slice - iMedian
            sliceMax = sliceMax - iMedan
            
        # Compute RMS
        vrange = vmax - vmin
        dvrange = vrange/10.
        maskRms1 = (v_slice >= vmin) & (v_slice <= (vmin+dvrange))
        maskRms2 = (v_slice >= (vmax-dvrange)) & (v_slice < vmax)
        rms1 = np.std(I_slice[maskRms1])
        rms2 = np.std(I_slice[maskRms2])
        rms = (rms1 + rms2)/2.
        individual_rms.append(rms)

        # Interpolate onto uniform velocity grid.
        # interp requires increaseing x axis.
        interp_I = np.interp(vel_grid, np.flip(v_slice), np.flip(I_slice))
        
#        print("I slice min, max: %.3f, %.3f" % (min(I_slice),max(I_slice)))
#        print("I inter min, max: %.3f, %.3f" % (min(interp_I),max(interp_I)))

        # Weighted sum
        summed += (w * interp_I)
        weightSum = weightSum + w
        # count plots to deal first update of offset
        nplot = nplot + 1

        if plot == 'both' or plot == 'stack':
            textOffset = plotOffset + (offset/3.)
            # Plot stacked spectrum
            plt.plot(v_slice, I_slice + plotOffset, lw=2)
            plt.text(vmin + 1., textOffset,
                     f"{nu_rest:.3f}", fontsize=10)
            plt.text(vmax - 3., textOffset,
                 labels[i], fontsize=10)

            # must give rest of plots some space.
            dOffset = max( offset, (max(I_slice)/2.))
            plotOffset = plotOffset + dOffset

        # end for all rest frequencies
        
    # now normalize weights for strongest line
    nSum = len(summed)
    print("Weight Sum: %.3f for %d samples" % (weightSum, nSum ))
    summed = summed * (1./ weightSum)
    
    if plot == 'both' or plot == 'sum':

        # now plot average
        plt.plot(vel_grid, summed + plotOffset, lw=3)
        # put text above plotted lines
        textOffset = plotOffset + (offset/3.)
        plt.text(vmax - 3., textOffset,
                 "Weighted Sum", fontsize=10)

        # plot the zero line for average
        plt.plot(vel_grid, (0.00001*summed) + plotOffset, 'g--', lw=.5)
        
    plt.xlabel("Velocity (km/s)", fontsize=14)
    if plot == 'sum':
        plt.ylabel("Average Intensity", fontsize=14)
    else:
        plt.ylabel("Intensity + offset", fontsize=14)
    # label for top of plot default
    if label == None:
        plt.title("Stack in Velocity %s with %s Spectra" %
                  (molecule, survey))
    else:
        plt.title(label)
    
    # Now fit a gaussian to the sum
    sumMax = summed.max()
    iMax = summed.argmax()
    velMax = vel_grid[iMax]
    rms = np.std(summed)

    # initial guess is at max location, with 1/20th toe velocity rangea
    initial_guess = [ sumMax, velMax, ((vmax-vmin)*.05)]
    try:
        popt, pcov = curve_fit(gaussian, vel_grid, summed, p0=initial_guess)
        # Extract the fitted parameters
        fit_A, fit_x0, fit_sigma = popt
        errors = np.sqrt(np.diag(pcov))
        print("Fitted parameters: Peak %10.6f +/- %10.6f" % (fit_A, errors[0]))
        print("                 : Velocity   %6.3f km/sec" % (fit_x0))
        print("                 : Line Width %6.3f km/sec" % (fit_sigma))
        # if velocity not specified and fit seems successful
        if velocity == None:
            vlabel = "Vel: %.3f$\pm$%.3f Width:%.3f$\pm$%.3f" % \
                (fit_x0, errors[1], fit_sigma, errors[2])
            velocity = fit_x0
        else:
            velocity = float(velocity)
            vlabel = "Vel: %.3f" % (velocity)
        print("Max Model Line   : %12.6f   %s (%.3f)" %
              (usedMaxFreq, usedMaxLabel, usedMaxWeight))
        plt.axvline(velocity, color='blue', linestyle="--", linewidth=1)
        plabel = "Peak: %.4f$\pm$%.4f" % (fit_A, errors[0])

        if plot == 'both' or plot == 'sum':
            plt.text(velocity+((vmax-vmin)*.05), plotOffset+sumMax,
                     vlabel, fontsize=10)
            plt.plot(vel_grid, gaussian(vel_grid, *popt)+plotOffset, 'r--')
            plt.text(velocity-((vmax-vmin)*.3), plotOffset+sumMax,
                     plabel, fontsize=10)
        else:
            plt.text(velocity+((vmax-vmin)*.05), plotOffset-dOffset+sliceMax,
                     vlabel, fontsize=10)
            # now do not double plot velocity
        velocity = None
        
    except RuntimeError as e:
        print(f"Error during fitting: {e}. No summed intensity detected.")

    if velocity != None:
        velocity = float(velocity)
        vlabel = "Vel: %.2f" % (velocity)
        plt.axvline(velocity, color='blue', linestyle="--", linewidth=1)
        plt.text(velocity+((vmax-vmin)*.01), plotOffset+sumMax,
                 vlabel, fontsize=10)

    # finally show result
    plt.tight_layout()
    plt.show()

    # create a file name with time and molecule without special characters
    datetime_iso = datetime.now().replace(microsecond=0).isoformat()
    moleculeNoDollar = molecule.replace("$","")
    moleculeNoLatex = moleculeNoDollar.replace("{","")
    moleculeNo = moleculeNoLatex.replace("}","")
    saveFile = "%s-%s.pdf" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    saveFile = "%s-%s.png" % (moleculeNo, datetime_iso)
    plt.savefig(saveFile)
    
    # Report S/N improvement
    if report_snr:
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
