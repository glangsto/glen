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

Using the spectral line stacking features are easy.   Below is an example for retrieving
and stacking GBT observations of Taurus Molecular Cloud and discover the radio Spectrum of the
cyanopolyene molecule HC5N.  The GOTHAM spectral line observations are downloaded on the first use
and placed in your ~/Downloads Directory

- git clone https://www.github.com/glangsto/glen
- cd glen
- ./glen -m pro/hc5n.pro
-  or
- ./glen --help
  
