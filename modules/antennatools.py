'''
  **Module with useful Antenna methods and definitions**
  
  :Author: S.F. Peik

'''
import matplotlib.pyplot as plt
from numpy import *

eta0 = 376.730


def Gain_to_AF(G,f):
    '''
    Parameters
    ----------
    G: float
        Gain in absolute
    f: float
        frequency in Hz
    
    Returns
    --------
    
    float:
       Antenna factor AF in absolut
       
    '''
    f = array(f)
    G = array(G)
    lam0 = 3e8/f
    AF = 2*pi/lam0*sqrt(2.4/G)
    return AF

def AF_to_Gain(AF,f):
    '''
    Parameters
    ----------
    AF: float
        AF in dB/m
    f: float
        frequency in Hz
    
    Returns
    --------
    
    float:
       Gain absolut
       
    '''
    f = array(f)
    AF = array(AF)
    lam0 = 3e8/f
    G = (9.73/lam0/10**(AF/20) )**2
    return G

def W_to_dBm(P):
    return 10*log10(P/1e-3)

def dBm_to_W(PdBm):
    return 10**(PdBm/10) * 1e-3

def dBm_to_dBuV(PdBm,Z0=50):
    return PdBm + 90 + 20*log10(sqrt(Z0))

def dBuV_to_dBm(PdBuV,Z0=50):
    return PdBuV - 90 - 20*log10(sqrt(Z0))

def W_to_dBuV(P, Z0=50):
    return 20*log10(sqrt(P*Z0)/1e-6)
    
def V_to_dBuV(V):
    return 20*log10(V/1e-6)
    
def lin_to_dB(x):
    return 10*log10(x)

def dB_to_lin(x):
    return 10**(x/10)


def friis(Pt: float,Gt,Gr,r,f)  -> float:
    '''
    Friis Transmission Formula in absolute values

    Parameters
    ----------
    Pt: Transmit power in W (float)
    Gt: Gain ransmit Antenne (float)
    Gr: Gain Receiving Antenna (float)
    r: Distance in m (float)
    f: frequency (float) 

    Returns
    -------

    float: 
       recevied power in W
    '''
    lam0 = 3e8/f
    return Pt*Gt*Gr*(lam0/4/pi/r)**2

def friis_dB(PtdBm: float,GtdBi,GrdBi,r,f)  -> float:
    '''
    Friis Transmission Formula in dB Values

    Parameters
    ----------
    Pt: Transmit power in dBm or dBW (float)
    Gt: Gain ransmit Antenne in dBi (float)
    Gr: Gain Receiving Antenna in dBi (float)
    r: Distance in m (float)
    f: frequency (float) 

    Returns
    -------

    float: 
       recevied power in W
    '''
    lam0 = 3e8/f
    LpdB = 20*np.log10(4*pi*r/lam0)
    return PtdBm+GtdBi+GrdBi-LpdB


from scipy.special import sici
def Si(x):
    return sici(x)[0]
def Ci(x): 
    return sici(x)[1]

def dipole_input_impedance(kl,a_to_l=0.0001):
    '''
    calculates the complex input impedance of a wire antenna
    kl: length in lambdas * 2pi    
    a_to_l: radius to diameter ratio of rod
    '''
    eta0 = 377.0
    game = 0.57721566499
    R = eta0/2/pi/(sin(kl/2))**2 * \
        (game +log(kl) - Ci(kl) 
         + 0.5 * sin(kl) * ( Si(2*kl) - 2*Si(kl))
         + 0.5 * cos(kl) * ( Ci(2*kl)-2*Ci(kl)+game+log(0.5*kl) ) 
        )
    
    X = eta0/2/pi/(sin(kl/2))**2 * \
        ( Si(kl) 
         + 0.5 * cos(kl) * ( -Si(2*kl) + 2*Si(kl))
         + 0.5 * sin(kl) * ( Ci(2*kl)-2*Ci(kl) + Ci(2*kl*(a_to_l)**2)) 
        )
    return R, X

def dipole_radiation_resistance(kl,a_to_l=0.0001):
    '''
    calculates the radiation resistance from input impedance of a wire antenna
    see Balanis 4-79
    kl: length in lambdas * 2pi    
    a_to_l: radius to diameter ratio of rod
    '''
    Rin, Xin = dipole_input_impedance(kl,a_to_l)
    Rr = Rin * (sin(kl/2))**2
    return Rr
    
    
    
def dipole_mutual_impedance_side_by_side(k,l,d):
    '''
    calculates the complex mutual impedance of a wire antenna
    see Balanis eqn 8.71 ff
    k: wave number 2pi/lam
    l: laenge   
    d: distance
    '''
    eta0 = 377.0
    u0 = k*d
    u1 = k*(sqrt(d**2+l**2)+l) 
    u2 = k*(sqrt(d**2+l**2)-l) 
    R21 =  eta0/4/pi*( 2*Ci(u0) - Ci(u1) - Ci(u2))
    X21 = -eta0/4/pi*( 2*Si(u0) - Si(u1) - Si(u2))
    
    return R21, X21


def AF_linear_array(theta, N, d = 0.5, delta_phase= 0 , AmpTaper = None):
    '''
    Array Factor of a linear array of N isotropic point sources, separated by distance d in lambda's  and
    all fed with constant magnitude and a constant phase increment psi
    --> Amplitude Taper not included yet!!
    :param theta: anglelist in degrees
    :param N: number of elements
    :param d: distance od elements in lambda's
    :param delta_phase: constant Phase difference between elements in degrees
    :param AmpTaper:
    :return:
    '''
    th = theta * pi/180
    lam = 1.0
    k = 2*pi/lam
    if AmpTaper is None:
        AmpTaper = ones(N)
    zeta = delta_phase * pi/180
    psi =  zeta + k*d * cos(th+pi/2) # hartnagel page 508
    #AF = abs(1/N * sin(N*psi/2) / sin(psi/2)) + 1e-6 # hartnagel page 508
    summed = 0
    for n in range(N):
        summed += exp(1j*n*psi) * AmpTaper[n]
    AF = abs( 1/N * summed )
    return AF

    
def polarPlotdB(theta,GdB,dBmin= -60, dBmax = 0, dBgrid = 10, minorGrid=True, angleSteps = 15, ax = None, **kwargs):

    '''
    :param theta: angle in degrees
    :param GdB: Gain or similar in dB values
    :param dBmin: minimum dB shown in center of plot
    :param dBmax: maximum dB shown on rim of plot
    :param dBgrid: grid of dB values
    :param minorGrid: show grid minor grid rings
    :param angleSteps: ticmarks at every angleStep degrees
    :return: fig,ax
    '''
    dBmin = int(dBmin)
    dBmax = int(dBmax)
    if ax is None:
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    thetarad = deg2rad(theta)
    ## The dBmin must be subtracted from all r values, as no neg. r values possible for polar plot

    ax.set_xlim(-pi, pi)
    ax.set_xticks(arange(-pi, pi, angleSteps * pi / 180))
    ## minor rings befor plot, such that it appears behind plot
    phi = arange(0, 2 * pi + pi / 100, pi / 100)
    for R in range(dBmin, dBmax, 2):
        rr = ones_like(phi) * R - dBmin
        ax.plot(phi, rr, "lightgray", lw=0.4)
    ## Plot now
    ax.plot(thetarad, GdB - dBmin, *kwargs.values())
    ## Set all adjustments
    ax.set_rmin(0.0)
    ax.set_rmax(-dBmin + dBmax)
    ax.set_rlabel_position(0)
    ## redefine for using nautical compass system ######
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ## Setting ticks ###################################
    ax.set_rticks(range(0, -dBmin + dBmax  + 1, dBgrid))  # Define the yticks
    labels = map(str, arange(dBmin, dBmax  + 1, dBgrid))
    dBlabels = [xxx + "dB" for xxx in labels]
    ax.set_yticklabels(dBlabels, fontsize='small')  # Change the labels
    ax.grid(True)
    try:
        return fig,ax
    except:
        return
        	    

if __name__ == '__main__':
    theta = arange(0,361,1.0)
    AF = AF_linear_array(theta,N=8, d=0.3)
    #G = sin(theta*pi/180)**2
    G = AF**2
    GdB = 10*log10(G+1e-10)
    fig, ax = polarPlotdB(theta,GdB,dBmin=-40,dBmax = 3, angleSteps=15, color="r")
    ## Select a sector of the graph
    ax.set_xlim(-pi/2,pi/2)
    ax.set_title("\n\nGain  $G(\\Theta)$")
    plt.show()
    plt.tight_layout()



