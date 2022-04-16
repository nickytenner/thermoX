from scipy import constants as cs
KB = cs.k              # Boltzmann constant kgm^2/s^2K
NA = cs.N_A            # Avogadro number
R  = NA*KB/4184.       # Universal gas constant kcal/molK
H  = cs.h              # kgm^2/s
AMU = 1.660538921e-27  # Atomic mass unit kg
BOHR = 5.292e-11       # Bohr radius m
C  = cs.c * 100        # speed of light in cm/s
search_strings = ["CARTESIAN COORDINATES (ANGSTROEM)",
                  "CARTESIAN COORDINATES (A.U.)",
                  "VIBRATIONAL FREQUENCIES",
                  "NORMAL MODES",
                  "FINAL SINGLE POINT ENERGY",
                  "Point Group:",
                  "The molecule is recognized as being linear",
                  "Final Gibbs free energy"]
