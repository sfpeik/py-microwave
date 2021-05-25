#
# Plot the response of a BFP420 Transistor
# uses the Touchstone file provided by Infineon
# S.F. Peik 2021
#


import mwave as mw
from numpy import *

f, S = mw.load_touchstone("bfp420.s2p")
fig, ax = mw.plotspar(f, S)
ax.set_ylim(-30, 30)
fig.show()

# Adding a Series Capacitor in the input
C1 = 1e-12
ABCD1 = ABCD4 = mw.ABCDseries(1/(1j * 2*pi*f * C1))  # add a series inductor
ABCD2 = mw.StoABCD(S, 50.0)
ABCD = mw.cascade([ABCD1, ABCD2])  # Cascade the five networks
Snew = mw.ABCDtoS(ABCD)
fig, ax = mw.plotspar(f, Snew)
ax.set_title("BFP420 with series Cap in Input line")
ax.set_ylim(-30, 30)
fig.show()
