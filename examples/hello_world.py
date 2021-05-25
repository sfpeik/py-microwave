# Minimal example to show that module mwave works
#
# S.F. Peik 2021
#

import mwave as mw

mw.hello()
print("Speed of Light is ", mw.c)
print("ABCD-matrix is: \n", mw.ABCDseries(100.0)) # ABCD MAtrix of a series 100 Ohm Resistor




