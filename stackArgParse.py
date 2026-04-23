#python
#HISTORY
#26Apr22 GIL Update minimum offset y axis offset
#26Apr21 GIL Allow multiple observation files
#26Feb23 GIL make ignoring confusing lines an option
#26Feb21 GIL add dv and baseline options
#26Feb20 GIL add usage and more options
#26Feb18 GIL initial version

import argparse
import textwrap

epilogText = textwrap.dedent('''
This program is intented for use in comparing molecule detections, based on 
FITS tables of multiple molecular line observations of spectral intensity 
versus frequency.   

The GOTHAM Spectral Pipeline data are available online at:
    https://greenbankobservatory.org/portal/gbt/gbt-legacy-archive/gotham-data/
References are there. Download the calibrated data.

The Turner and Langston Q band Survey of TMC-1 and searches for '$HC_{13}N$' are
in preparation, once the spectral comparision is complete.   These
data will also be downloadable.

See Langston and Turner (2007) for a first example of successful molecular 
line stacking to detect previously un-detected molecular species. 

Article: "Detection of C Isotopomers of the Molecule HC_7N" 
(Langston, G., & Turner, B., 2007, The Astrophysical Journal, 658, 455).

Example:
 ./glen -i pro/hc5n.pro --molecule '$HC_{5}N$' --plot 'sum'

Line summing code crafted by Glen Langstion, 2026 February 23
'''
                             )
# parse the extensive argument list for stacking spectra
def stackArgParse():
    parser = argparse.ArgumentParser(

        description="Process molecular spectra from one or more data files.   Spectral lines are summed in velocity space, using modeled weighted averages.",
        epilog=epilogText,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )


    # Optional distance to source, for calculation of molecule mass
    parser.add_argument(
        "-b", "--baseline",
        action='store_true',
        help="Optionally subtract median of data as baseline constant"
    )

    parser.add_argument(
        "-c", "--cologne",
        type=str,
        help="Cologne Molecular Spectral intensity model file (not yet implemented)"
    )

    # One or more optional data files
    parser.add_argument(
        "-d", "--data-files",
        nargs="+",          # one or many
        type=str,
        default="tmc-tlq.fits gotham_drv.fits",
        help="Input data files: will look in data directory, ie tmc-tql.fits"
    )

    # Optional string argument
    parser.add_argument(
        "--ignore",
        type=str,
        help="Ignore observations of confusing lines not in center of Sum Range"
    )

    # Optional string argument
    parser.add_argument(
        "-f", "--frequency",
        type=str,
        help="Rest frequency of Molecule (MHz)"
    )

    # Optional string argument
    parser.add_argument(
        "-g", "--gauss",
        action='store_true',
        help="Fit Gaussians to Lines"
    )

    # Optional distance to source, for calculation of molecule mass
    parser.add_argument(
        "-i", "--intensity",
        type=str,
        help="File containing the line intensities model in IDL format (historical)"
    )

    # Optional distance to source, for calculation of molecule mass
    parser.add_argument(
        "-l", "--light-years",
        type=str,
        help="Distance to source, in light years, used to calculate total molecular mass"
    )

    # Optional string argument
    parser.add_argument(
        "-m", "--molecule",
        type=str,
        help="Name of the molecule (e.g., HC7N or '$HC_7N$')"
    )

    # Optional string argument
    parser.add_argument(
        "-n", "--normalize",
        type=str,
        help="Type of intensity normalization: none, peak or rms"
    )

    # Optional string argument
    parser.add_argument(
        "-o", "--offset",
        type=str,
        default = "0.5",
        help="Intensity offset between plots (K)"
    )

    # Optional string argument
    parser.add_argument(
        "-p", "--plot",
        type=str,
        default='both',
        help="Type of plot: one of 'stack', 'sum', 'freq' or 'both'. Default 'both'"
    )

    # Optionally report SNR of fit
    parser.add_argument(
        "--report-snr",
        action='store_true',
        help="Print summary of Fit to stacked line."
    )

    # Optional Survey Name for label
    parser.add_argument(
        "-s", "--survey",
        type=str,
        default="TMC-TLQ GOTHAM",
        help="Spectral Observation Survey Name."
    )

    # Optional physical temperature of molecule model
    parser.add_argument(
        "-t", "--temperature",
        type=str,
        help="Model Temperature used for calculations of the spectral intensities"
    )

    # Optional string argument
    parser.add_argument(
        "-v", "--velocity",
        type=str,
        help="Relative Velocity of object (km/sec, LSRK)"
    )

    # Optional string argument
    parser.add_argument(
        "-vmin",
        type=str,
        default="1.",
        help="Minimum Velocity for summed spectrum"
    )

    # Optional string argument
    parser.add_argument(
        "-vmax",
        type=str,
        default="11.",
        help="Maximum Velocity for summed spectrum"
    )

    parser.add_argument(
        "-dv", 
        type=str,
        default="0.05",
        help="Velocity increment (km/sec) for gridded summed spectra"
    )

    # Optional absolute weight of strongest predicted line in line list
    parser.add_argument(
        "-w", "--weight",
        type=str,
        help="Absolute weight of strongest line in line list.  Used for abundance calculation"
    )

    # Example of another optional string
    parser.add_argument(
        "--title",
        type=str,
        help="Optional title for top of output plots"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = stackArgParse()
    print("Molecule:", args.molecule)
    print("Data files:", args.data_files)
    print("Label:", args.title)
    print("Velocity:", args.velocity)
    if args.frequency != None:
        print("Frequency:", args.frequency)
