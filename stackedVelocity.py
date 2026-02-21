#python
#HISTORY
#26Feb15 GIL initial stacked velocity plot

import numpy as np
import matplotlib.pyplot as plt

# Speed of light in km/s
C = 299792.458

def freq_to_velocity(nu_obs, nu_rest):
    """
    Convert observed frequencies to velocities using the radio definition:
        v = c * (nu_rest - nu_obs) / nu_rest
    """
    return C * (nu_rest - nu_obs) / nu_rest


def stacked_velocity_plots(freq, intensity, rest_freqs, labels,
                           vmin=0.0, vmax=12.0, offset=0.5):
    """
    Create stacked velocity plots for a set of molecular line rest frequencies.

    Parameters
    ----------
    freq : array
        Observed frequency array (same units as rest_freqs)
    intensity : array
        Observed intensity array
    rest_freqs : list or array
        Rest frequencies for the molecular transitions
    labels : list
        Text strings for labeling plots
    vmin, vmax : float
        Velocity window to extract (km/s)
    offset : float
        Vertical spacing between stacked spectra
    """

    plt.figure(figsize=(8, 10))

    nplot = 0
    lastMax = 0.
    plotOffset = 0.
    for i, nu_rest in enumerate(rest_freqs):
        # Convert entire spectrum to velocity for this line
        vel = freq_to_velocity(freq, nu_rest)

        # Select velocity window
        mask = (vel >= vmin) & (vel <= vmax)
        v_slice = vel[mask]
        I_slice = intensity[mask]

        if len( v_slice) < 10:
            print("Not enough samples for frequency %10.3f" % \
                  (nu_rest))
            continue

        # accumulate offsets 
        plotOffset = plotOffset + (lastMax*.5)
        # Apply vertical offset for stacking
        I_offset = I_slice + plotOffset
        nplot = nplot + 1
        # add an additional bit of an offset for next plot
        lastMax = max(I_slice)

        # Plot the stacked spectrum
        plt.plot(v_slice, I_offset, lw=2)

        textOffset = plotOffset + (offset*.333)
        # Label each stacked spectrum
        plt.text(vmin + 0.2, textOffset,
                 f"{nu_rest:.3f}", fontsize=10)
        plt.text(vmax - 3 , textOffset,
                 labels[i], fontsize=10)

        # prepare for next plot
        plotOffset = plotOffset + max( lastMax*.5, offset)
        
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Intensity + offset")
    parts = labels[0].split(" ")
    molecule = parts[0]
    plt.title("Stacked %s Spectra in Velocity Space" % (molecule))
    plt.xlim(vmin, vmax)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------
# Example usage (replace with your real data)
# ---------------------------------------------------------
if __name__ == "__main__":
    # Example synthetic data
    freq = np.linspace(45e9, 46e9, 50000)       # 45–46 GHz
    lables = [ "line 1", "line 2", "line 3"]
    intensity = np.random.normal(0, 0.02, len(freq))

    # Example rest frequencies for a molecule
    rest_freqs = [
        45.123e9,
        45.456e9,
        45.789e9
    ]

    stacked_velocity_plots(freq, intensity, rest_freqs, labels,
                           vmin=0, vmax=12, offset=0.3)
