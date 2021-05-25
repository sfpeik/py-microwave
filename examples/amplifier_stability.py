#
# Check an amplifier or transistor for stability
#
#
import mwave as mw
from numpy import *
from matplotlib.pyplot import show

# Load data from bfp420 transistor
flist, Slist = mw.load_touchstone("bfp420.s2p")
# find entry at frequency
fwanted = 300e6
idx = abs(flist-fwanted).argmin()
f, S = flist[idx], Slist[idx]
print(f, '\n', S)

## Calculate stability data
Cs,Rs,Cl,Rl,mu1, mu2, fig, ax = mw.AmpStabilityCircle(S, plotit=True)
print("mu1=",mu1)
print("mu2=",mu2)
print(f"Source Stab. Circle at: C={mw.magphase_str(Cs)}  R={Rs:5.3f}")
print(f"Load   Stab. Circle at: C={mw.magphase_str(Cl)}  R={Rl:5.3f}")
show()

