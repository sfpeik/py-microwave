#
# Plot the response of a Low-Pass Filter S-Parameters
#
# S.F. Peik 2021
#
# We build up a common lowpass in ladder structure with
#
#  o--o--(((--o--(((--o--o
#     |       |       |
#    ===     ===     ===
#     |       |       |
#    ---     ---     ___
#

import mwave as mw
from numpy import *

### Set values according to filter theory ######
### Using G-values for Chebychev filter with 0.5 dB Ripple
### and corner freq fc
Z0 = 50.0 ## ref. Port Impedanace
fc = 2e9 # corner frequency
wc = 2*pi*fc
C1 = 1.7058 / wc / Z0
L2 = 1.2296/ wc * Z0
C3 = 2.5409/ wc / Z0
L4 = L2
C5 = C1

mw.hello()
f = arange(1e6,4e9,1e6)
w = 2 * pi * f
ABCD1 = mw.ABCDshunt(1j * w * C1)  # add a shunt capacitor
ABCD2 = mw.ABCDseries(1j * w * L2)  # ad a series inductor
ABCD3 = mw.ABCDshunt(1j * w * C3)  # add a shunt capacitor
ABCD4 = mw.ABCDseries(1j * w * L4)  # ad a series inductor
ABCD5 = mw.ABCDshunt(1j * w * C5)  # add a shunt capacitor
ABCD = mw.cascade([ABCD1, ABCD2, ABCD3, ABCD4, ABCD5])  # Cascade the five networks
S = mw.ABCDtoS(ABCD,Z0)
fig, ax = mw.plotspar(f,S)
ax.set_ylim(-30,2)
fig.show()



