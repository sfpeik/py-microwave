import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import mwave as mw

# Create a simple Smtih chart
# this will be superimposed over the generated plot.
#
def plot_smith(ax):
    '''
    Simple Smith Chart
    '''
    for Z_real  in (0,0.33,1,3):
        x = np.zeros(100)
        y = np.zeros(100)
        for ii in range(100):
            Zc = Z_real+1j*np.exp(ii/10.-5)
            Gam = (Zc-1)/(Zc+1)
            x[ii] = np.angle(Gam)
            y[ii] = np.abs(Gam)
        x[0] = np.pi
        x[-1] = 0
        ax.plot(x,y,'k-')
        ax.plot(-x,y,'k-')
    for Z_imag  in (0,0.33,1.,3.):
        x = np.zeros(100)
        y = np.zeros(100)
        for ii in range(100):
            Zc = np.exp(ii/10-5) + 1j*Z_imag
            Gam = (Zc-1)/(Zc+1)
            x[ii] = np.angle(Gam)
            y[ii] = np.abs(Gam)
        y[0] = 1
        x[-1] = np.pi
        ax.plot(x, y,'k-',zorder=2)
        ax.plot(-x, y,'k-')
        ax.set_yticks([])
        ax.grid(False)

'''
The following plot shows the variation of the transducer gain 
when chosing the source impedance of the amplifier. 
Clearly seen, the gain is maximised when choosing Γs to be conj(S11)
'''
fig,ax = plt.subplots(figsize=(10,10),
subplot_kw=dict(projection='polar'))
mag = np.linspace(0, 1, 100)
phase = np.linspace(0, 2 * np.pi, 100)
r, th = np.meshgrid(mag, phase)
Gams = r*np.exp(1j*th)
S = np.matrix([[0.3+0.6j,0.01],[10,0.4j]])
Gt = mw.transducerGain(S,Gams,0.0)
myplot = ax.pcolormesh(th, r, Gt,
       cmap='YlOrRd',shading='gouraud',
       vmin=15., vmax=24)
ax.set_xticklabels([])
fig.colorbar(myplot)
# Plot a point at the conj. S11
Sconj = np.conj(S[0,0])
ax.plot(np.angle(Sconj), np.abs(Sconj), "w*", ms=20)
plot_smith(ax)
fig.savefig('smith_trans_gain.png')
fig.show()