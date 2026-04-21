#/bin/python
#Fit gaussians to a spectrum
#HISTORY
#26Apr11 GIL fix three_gaussians to deal with single float versus array
#26Apr10 GIL/Copilot initial version

import numpy as np
from scipy.optimize import curve_fit

def gaussian(x, amp, cen, wid):
    """Single Gaussian."""
    return amp * np.exp(-(x - cen)**2 / (2 * wid**2))

def three_gaussians(x, *params):
    """Sum of up to three Gaussians."""
    g = np.zeros_like(x)
    if not isinstance( params[0], np.float64):
        params = params[0]
# diagnostics to figure out pointer problems
#    print("Params: %d" % (len(params)))
#    print(params)
    
    for i in range(0, len(params), 3):
        amp, cen, wid = params[i:i+3]
#        print("%3d: %.2f %.2f %.2f" % (i, amp, cen, wid))
        g += gaussian(x, amp, cen, wid)
    return g


def fit_single_gaussian(freq, inten, guess=None):
    """Fit one Gaussian to the data."""
    if guess is None:
        # amplitude guess = max intensity
        amp0 = np.max(inten)
        cen0 = freq[np.argmax(inten)]
        wid0 = (freq[-1] - freq[0]) / 20.0
        guess = [amp0, cen0, wid0]

    try:
        popt, pcov = curve_fit(gaussian, freq, inten, p0=guess)
    except RuntimeError:
        print("Single Gaussian fit failed.")
        popt = [0, 0, 1]
    return popt


def iterative_three_gaussian_fit(freq, inten):
    """
    Iteratively:
    1. Fit Gaussian to full spectrum
    2. Subtract it
    3. Find next peak
    4. Fit second Gaussian
    5. Repeat for third
    6. Fit all three simultaneously
    """
    residual = inten.copy()
    fits = []

    print("\n=== Iterative Gaussian Fitting ===")

    for n in range(3):
        print(f"\n--- Fitting Gaussian #{n+1} ---")

        # Find peak in residual
        peak_idx = np.argmax(residual)
        cen_guess = freq[peak_idx]
        amp_guess = residual[peak_idx]
        wid_guess = (freq[-1] - freq[0]) / 30.0

        guess = [amp_guess, cen_guess, wid_guess]
        print(f"Initial guess: amp={amp_guess:.3f}, cen={cen_guess:.6f}, wid={wid_guess:.6f}")

        popt = fit_single_gaussian(freq, residual, guess)
        fits.append(popt)

        print(f"Fit #{n+1}: amp={popt[0]:.3f}, cen={popt[1]:.6f}, wid={popt[2]:.6f}")

        # Subtract fitted Gaussian
        residual = residual - gaussian(freq, *popt)

    print("\n=== Final 3-Gaussian Combined Fit ===")

    # Flatten initial guesses for combined fit
    p0 = np.array(fits).flatten()

    try:
        popt_all, pcov_all = curve_fit(three_gaussians, freq, inten, p0=p0)
    except RuntimeError:
        print("Final 3-Gaussian fit failed.")
        return fits, None

    print("\nFinal combined fit parameters:")
    for i in range(3):
        amp, cen, wid = popt_all[3*i:3*i+3]
        print(f"Gaussian {i+1}: amp={amp:.3f}, cen={cen:.6f}, wid={wid:.6f}")

    return fits, popt_all


# Example usage:
if __name__ == "__main__":
    # Example synthetic data
    freq = np.linspace(10.0, 11.0, 2000)
    true_params = [1.0, 10.30, 0.01,
                   0.7, 10.55, 0.015,
                   0.4, 10.80, 0.02]
    inten = three_gaussians(freq, true_params)
    inten += 0.05 * np.random.randn(len(freq))  # add noise

    fits, final = iterative_three_gaussian_fit(freq, inten)
