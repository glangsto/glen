# GLEN — GBT Line Emission Navigator

**Navigating faint molecular signals through stacked spectral insight.**

GLEN is a high‑performance spectral‑line stacking and molecule‑discovery engine designed for the Green Bank Telescope (GBT). It converts measured spectra into velocity space, aligns multiple molecular transitions, applies optional weighting and normalization, and produces both stacked and summed profiles with optimized signal‑to‑noise.

GLEN is built for precision molecular searches, enabling the detection of faint, blended, or marginally detected lines through coherent velocity‑space integration.

---

## Features

- **Velocity‑space conversion** for each molecular transition  
- **Stacked spectral visualization** with customizable offsets  
- **Weighted stacking** for line‑strength or S/N‑based weighting  
- **Normalization modes**: peak, area, RMS, or none  
- **Baseline subtraction**  
- **Gaussian fitting** of the summed spectrum  
- **S/N improvement diagnostics**  
- **Uniform velocity‑grid interpolation** to avoid sampling gaps  
- **Clean, publication‑quality plots**  

---

## Why GLEN?

Modern molecular searches often rely on stacking multiple weak transitions to reveal faint emission. GLEN provides a robust, flexible, and scientifically rigorous framework for:

- detecting new molecular species  
- improving S/N through coherent stacking  
- comparing transitions across a molecule’s ladder  
- producing reproducible, publication‑ready figures  

---

## Quick Start

```python
from glen import stackSumVelocity

vel_grid, summed, fit_params = stackSumVelocity(
    freq, intensity, rest_freqs,
    vmin=0, vmax=12,
    dv=0.05,
    offset=0.3,
    weights=[1, 0.9, 0.8, 0.7],
    normalize="peak",
    baseline_subtract=True,
    fit_gaussian=True,
    report_snr=True
)
