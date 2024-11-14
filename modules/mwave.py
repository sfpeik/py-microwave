##### use by:
## import sys
## sys.path.append('/home/speik/soridat/pythonlib')
## import mwave as mw
##
### Uses Napoleon extension for Docstring generation, see https://sphinxcontrib-napoleon.readthedocs.io/en/latest/

'''
  **Module with useful Microwave methods and definitions**
  
  :Author: S.F. Peik

'''


__license__ = "Soeren Peik"
__revision__ = " 2021-April "
__docformat__ = 'reStructuredText'
__version__ = "1.0.1"


from numpy import  array,sqrt,pi,log,matrix, conj, angle, zeros, exp, abs, ndim, log10, arange, around, \
    shape, ones, tan, isnan, nan, cosh, sinh, atleast_1d, transpose, squeeze, zeros_like, ones_like, \
    broadcast_to, identity, matrix, real, imag, tanh, interp, set_printoptions, asarray, argmin, concatenate, newaxis
import matplotlib.pyplot as plt
import sys

try:
    import smith as smi
except:
    print("Smith Module not found, doing without")

from scipy.optimize import fsolve, brentq, brenth, brent

#### some constants

#: Permeability of free space
mu0  = 4*pi*1e-7   
#: velocity of light
c    = 299792458.0  
#: Permittivity of free space
eps0 = 1/mu0/c/c    
#: Boltzman Constant
k    = 1.3806488e-23 
#: Free Space impedance
eta0 = sqrt(mu0/eps0) 

def hello():
    print("Here is py-microwave Version:", __version__)
########################################################################
def coth(x):
    return 1/tanh(x)

def find_nearest(arr, value):
    '''
    find the index of a nearest value in an array
    '''
    arr = asarray(arr)
    idx = (abs(arr - value)).argmin()
    return idx

########################################################################    
def lineinputimpedance(Z0,Zl,betal):
    r'''Calculates input impedance of a terminated line
    
    Parameters
    ----------
    Z0 : complex 
        Line impedance in Ohm (type complex)
    Zl : complex
        Load impedance in Ohm (type complex)
    betal : float
        Electrical length :math:`\beta l` in radians (type float)
    
    Returns 
    -------
    complex
        Input impedance in Ohm 
    
    Examples
    --------
    >>> Zin = lineinputimpedance(50,100,3.14/2)
    >>> around(Zin)
    (25-0j)
    '''
    
    Zin = Z0 * (Zl+1j*Z0*tan(betal)) / (Z0+1j*Zl*tan(betal))
    return Zin

###########################################################################
def msimpedance(w,h,er):
  r'''
    Calculates microstrip line impedance :math:`Z_0` and :math:`\epsilon_{eff}` from Wheeler formula 
    
    Parameters
    ----------
    
    w : float
        width of Microstrip-line 
    h : float
        height of substrate 
    er : float
        Epsilon relative of substrate
    
    Returns
    -------
    
    tuple
        Z0, eps_eff 
     
    Examples
    --------    
    
    >>> msimpedance(1.8,0.83,3.55)
    (51.129452787061204, 2.773818757727919)
  '''

  if w ==0: w=1e-12
  eta0=377
  if(w/h<=1):
      F=1/sqrt(1+12*h/w)+0.004*(1-w/h)**2
  else:
      F=1/sqrt(1+12*h/w)
  e_eff=0.5*(er+1+(er-1)*F)

  if(w/h<=1): 
      Z0=eta0/sqrt(e_eff)*1/2/pi*log(8*h/w+0.25*w/h)
  else:
      Z0=eta0/sqrt(e_eff)*1/(w/h+2.46-0.49*h/w+(1-h/w)**6)
  return Z0,e_eff

###########################################################################
def msdimension(Z0wanted,elen,f,h,epsr):

    r'''
    Calculates microstrip line dimensions from the impedance, elec. length
    
    Note
    ----
    This function uses an optimizer to find the dimensions
    
    Parameters
    -----------
    
    Z0wanted : float
        Impedance of MS-line 
    elen : float
        elec length of line in :math:`\lambda` 
    f : float
        frequency 
    h : float
        height of substrate in mm 
    epsr : float
        Epsilon relative of substrate 
        
    Returns  
    --------
    
    tuple
        w,l, epseff 
    
    
    Examples
    --------
    Here is an Example
    
    >>> msdimension(50,1/4.,2e9,0.83,3.55)
    (1.867, 0.02248463090701959, 2.7815765073092837)
    >>> msdimension(2,1/4.,2e9,0.83,3.55)
    Traceback (most recent call last):
        ...  
    ValueError: could not find solution in msdimension for 2.000000 Ohms
    '''

    lam0=3e8/f
    imp= lambda w:  msimpedance(w,h,epsr)[0]-Z0wanted
    try:
        result = brentq(imp,0.002*h,20.0*h, xtol=1e-15) 
    except:
        raise ValueError('could not find solution in msdimension for %f Ohms' % (Z0wanted))
    w=round(result,8)
    Z0,epseff=msimpedance(w,h,epsr)
    lamms=lam0/sqrt(epseff)
    l = elen * lamms
    return (w,l,epseff)
 
########################################################################   
def coupledlineCaps(w,h,s,epsr):
    '''
    used for the calculation of msCoupledLineImp
    '''
    
    Zc,epseff=msimpedance(w,h,epsr)
    #print(Zc,epseff)
    Cp   = eps0*epsr*w/h             # 4.24
    Cf   = 1/2. * (sqrt(epseff)/(c*Zc) - Cp)  # 4.26
    A =exp(-0.1*exp(2.33-2.53*w/h))     # 4.xx
    Cfpr = Cf / ( 1+A*(h/s) *tanh(8*s/h))   # 4.26

    Cgd = eps0*epsr/pi * log( coth(pi/4*s/h) ) + 0.65* Cf * (0.02*sqrt(epsr)/(s/h) + 1 - 1/epsr**2)

    k = (s/h) / (s/h + 2*w/h)
    
    kpr = sqrt(1-k**2)

    if     0<k**2 and k**2<0.5:
        K_to_K = 1/pi* log(2* (1+sqrt(kpr)) / (1-sqrt(kpr)) ) 
    elif 0.5<k**2 and k**2<1:
        K_to_K = pi / log(2* (1+sqrt(k)) / (1-sqrt(k)) ) 
    else:
        print("Outside Range")
        exit(1)

    Cga = eps0*K_to_K

    Ce = Cp + Cf + Cfpr
    Co = Cp + Cf + Cgd + Cga
    return Ce,Co

########################################################################

def msCoupledLineImp(w,h,s,epsr):
    '''
    Calculates microstrip even and odd impedances of a coupled microstrip line 
    
    Parameters
    ----------
    
    w : float
        width of Microstrip-line 
    h : float
        height of substrate 
    s : float
        gap between lines         
    er : float
        Epsilon relative of substrate
    
    Returns
    -------
    
    (float, float)
        even and odd impedance of coupled lines
     
    Examples
    --------    
    
    >>> msCoupledLineImp(1.5,0.5,0.2,2.56)
    (55.6657732006762, 36.95927306032294)
    '''
    Ce,Co = coupledlineCaps(w,h,s,epsr)
    Cae,Cao = coupledlineCaps(w,h,s,1.0)
    Ze = (c * sqrt(Cae*Ce))**-1  # 4.29
    Zo = (c * sqrt(Cao*Co))**-1  # 4.30
    return Ze,Zo

########################################################################

def msVia(h,D):
    '''
    Calculates inductance of a MS Via

    Parameters
    ----------

    h : float
        height of substrate
    D : float
        Diameter of via

    Returns
    -------

    float
        Inductance of via in Henry

    Examples
    --------

    >>> msVia(1.5e-3,0.9e-3)
    2.408551487285845e-10
    '''

    r = D/2
    L = mu0/2/pi * ( h*log( (h + sqrt(r*r+h*h)) / r )  + 3/2 * (r-sqrt(r*r+h*h)) )
    return L


### ABCD Matrix for Series Element #####################################
def ABCDseries(Z):
    r'''
    Creates an ABCD matrix (2x2 Array) for a series impedance element
    
    :math:`\left[ \begin{matrix} 1 & Z \\ 0 & 1 \end{matrix} \right]`
    
    Supports array type Z

    Examples
    --------
    
    >>> ABCDseries(100.0)
    array([[  1., 100.],
           [  0.,   1.]])
    >>> Z = arange(10,40,10)
    >>> ABCDseries(Z)
    array([[[ 1, 10],
            [ 0,  1]],
    <BLANKLINE>
           [[ 1, 20],
            [ 0,  1]],
    <BLANKLINE>
           [[ 1, 30],
            [ 0,  1]]])
    '''
    Z = atleast_1d(Z)
    x = array([[ones_like(Z),Z],[zeros_like(Z),ones_like(Z)]])
    x = transpose(x,(2,0,1))
    return squeeze(x)

### ABCD Matrix for Shunt Element #########################
def ABCDshunt(Y):
    Y = atleast_1d(Y)
    x = array([[ones_like(Y),zeros_like(Y)],[Y,ones_like(Y)]])
    x = transpose(x,(2,0,1))
    return squeeze(x)

### ABCD Matrix for Shunt Element #########################
def ABCDinverse(ABCD):
    x = []
    for AA in ABCD:
        x.append( matrix(AA)**-1 )
    return squeeze(x)

### ABCD Matrix for Line Element #########################
def ABCDline(beta, length,Z0,alpha=0.0):
    r'''
    Creates an ABCD matrix (2x2 Array) for an inserted line  element
       
    Parameters
    ----------
    
    beta : float
        phase constant at operating freq 
    l : float
        length in meter 
    Z0 : float
        line impedance       
    alpha : float
        loss in Np
    
    Returns
    -------
    
    2x2 matrix
        ABCD matrix at freq point beta
     
    
    Supports array type Z

    Examples
    --------
    
    >>> f = 2e9
    >>> lam = c/f
    >>> beta = 2*pi/lam
    >>> ABCDline(beta,0.1,50.0)
    array([[-0.49748657-0.00000000e+00j,  0.        -4.33735840e+01j],
           [ 0.        -1.73494336e-02j, -0.49748657-0.00000000e+00j]])

    '''
    beta = atleast_1d(beta)
    alpha = atleast_1d(alpha)
    gammal = (alpha+1j*beta)*length
    x = array([[cosh(gammal), Z0*sinh(gammal)],[ 1./Z0*sinh(gammal), cosh(gammal)]])
    x = transpose(x,(2,0,1))
    return squeeze(x)

############################################################
def cascade(ABCDlist):
    '''
    Cascades a list of 2x2 matrices, usually ABCD matrices
    
    Examples
    --------
    >>> ABCDlist = []
    >>> ABCDlist.append( ABCDseries(100) )
    >>> ABCDlist.append( ABCDshunt(20) )
    >>> ABCDlist.append( ABCDseries(50j) )
    >>> cascade(ABCDlist)
    array([[2.001e+03     +0.j, 1.000e+02+100050.j],
           [2.000e+01     +0.j, 1.000e+00  +1000.j]])
    
    '''

    try:
      if ndim(ABCDlist) == 3:
          fpoints = 0
          blocks,two,two = shape(ABCDlist)
          A = identity(2)
      else:
          blocks,fpoints,two,two = shape(ABCDlist)
          A = broadcast_to(identity(2),(fpoints,2,2))
      for ABCD in ABCDlist:
          A = A @ ABCD
      #print("Cascaded ",blocks," ABCD-Blocks")
      return A
    except: 
        print("in cascade: array not of shape n_blocks,freqpoints,2,2")



############################################################
def ABCDJInverter(J):
    r'''
    Create the ABCD Matrix of a J inverter with J as input
    
    :math:`\left[ \begin{matrix} 0 & j/J \\ j\cdot J & 0 \end{matrix} \right]`
    
    Examples
    --------
    
    >>> ABCDJInverter(5)
    array([[0.+0.j , 0.+0.2j],
           [0.+5.j , 0.+0.j ]])
    '''
    J = atleast_1d(J)
    x = array([[zeros_like(J),-1/(1j*J)],[1j*J,zeros_like(J)]])
    x = transpose(x,(2,0,1))
    return squeeze(x)

########################################################################
def ABCDKInverter(K):
    r'''
    Create the ABCD Matrix of a K inverter with K as input
    
    :math:`\left[ \begin{matrix} 0 & j\cdot K \\ j / K & 0 \end{matrix} \right]`
    
    Examples
    --------
    
    >>> ABCDKInverter(0.2)
    array([[0.+0.j , 0.+0.2j],
           [0.+5.j , 0.+0.j ]])
    '''
    K = atleast_1d(K)
    x = array([[zeros_like(K),1j*K],[-1/(1j*K),zeros_like(K)]])
    x = transpose(x,(2,0,1))
    return squeeze(x)


########################################################################
def ABCDtoS(ABCD,Z0=50):
    if shape(ABCD) == (2,2):
        ABCD = array([ABCD])
    A=ABCD[:,0,0]
    B=ABCD[:,0,1]
    C=ABCD[:,1,0]
    D=ABCD[:,1,1]
    S11=(A+B/Z0-C*Z0-D)/(A+B/Z0+C*Z0+D)
    S12=2*(A*D-B*C)/(A+B/Z0+C*Z0+D)
    S21=2/(A+B/Z0+C*Z0+D)
    S22=(-A+B/Z0-C*Z0+D)/(A+B/Z0+C*Z0+D)
    S=array([[S11,S12],[S21,S22]])
    S=transpose(S,(2,0,1))
    return squeeze(S)

############################################################
def Sinterpolate(fnewlist, flist, S):
    '''
    
    >>> flist = array([1.0e9,1.5e9,2.0e9])
    >>> fnewlist = arange(1.0e9,2.0e9,0.1e9)
    >>> ABCD = ABCDseries(2j*pi*flist*10e-9)  # create an array of ABCD matrices for a 10nH inductror
    >>> S = ABCDtoS(ABCD,50.0)
    >>> Sinterpolate(fnewlist, flist, S)
    array([[[0.2830432 +0.45047724j, 0.7169568 -0.45047724j],
            [0.7169568 -0.45047724j, 0.2830432 +0.45047724j]],
    <BLANKLINE>
           [[0.32051719+0.46020656j, 0.67948281-0.46020656j],
            [0.67948281-0.46020656j, 0.32051719+0.46020656j]],
    <BLANKLINE>
           [[0.35799118+0.46993589j, 0.64200882-0.46993589j],
            [0.64200882-0.46993589j, 0.35799118+0.46993589j]],
    <BLANKLINE>
           [[0.39546517+0.47966521j, 0.60453483-0.47966521j],
            [0.60453483-0.47966521j, 0.39546517+0.47966521j]],
    <BLANKLINE>
           [[0.43293916+0.48939453j, 0.56706084-0.48939453j],
            [0.56706084-0.48939453j, 0.43293916+0.48939453j]],
    <BLANKLINE>
           [[0.47041315+0.49912385j, 0.52958685-0.49912385j],
            [0.52958685-0.49912385j, 0.47041315+0.49912385j]],
    <BLANKLINE>
           [[0.49878519+0.49674541j, 0.50121481-0.49674541j],
            [0.50121481-0.49674541j, 0.49878519+0.49674541j]],
    <BLANKLINE>
           [[0.52715723+0.49436697j, 0.47284277-0.49436697j],
            [0.47284277-0.49436697j, 0.52715723+0.49436697j]],
    <BLANKLINE>
           [[0.55552928+0.49198854j, 0.44447072-0.49198854j],
            [0.44447072-0.49198854j, 0.55552928+0.49198854j]],
    <BLANKLINE>
           [[0.58390132+0.4896101j , 0.41609868-0.4896101j ],
            [0.41609868-0.4896101j , 0.58390132+0.4896101j ]]])
    '''   
    S11=S[:,0,0]
    S12=S[:,0,1]
    S21=S[:,1,0]
    S22=S[:,1,1]
    
    S11i = interp(fnewlist, flist, S11)
    S12i = interp(fnewlist, flist, S12)
    S21i = interp(fnewlist, flist, S21)
    S22i = interp(fnewlist, flist, S22)
    Si=array([[S11i,S12i],[S21i,S22i]])
    Si=transpose(Si,(2,0,1))
    return squeeze(Si)
    

############################################################
def Svalue(f,flist,S):
    '''
    interpolates the S-matrix for a given freq piont
    
    f: freq. point to be interpolated
    flist: array of freq points (float)
    S: array of (2,2) matrices same size as flist
    return:
    S-matrix for point f
    '''
    S11=S[:,0,0];S12=S[:,0,1];S21=S[:,1,0];S22=S[:,1,1]
    if len(flist) != len(S11):
        raise ValueError("Diffrent Length of f and S")
    Si = zeros((2,2),dtype=complex)
    Si[0,0] = interp(f,flist,S11)
    Si[0,1] = interp(f,flist,S12)
    Si[1,0] = interp(f,flist,S21)
    Si[1,1] = interp(f,flist,S22)
    return Si

############################################################
def StoABCD(S,Z0=50):
    if shape(S) == (2,2):
        S = array([S])
   
    S11=S[:,0,0]
    S12=S[:,0,1]
    S21=S[:,1,0]
    S22=S[:,1,1]
    A =      ((1+S11)*(1-S22)+S12*S21) / 2/S21
    B = Z0*  ((1+S11)*(1+S22)-S12*S21) / 2/S21
    C = 1/Z0*((1-S11)*(1-S22)-S12*S21) / 2/S21
    D =      ((1-S11)*(1+S22)+S12*S21) / 2/S21
    ABCD=array([[A,B],[C,D]])
    ABCD=transpose(ABCD,(2,0,1))
    return squeeze(ABCD)

#############################################################
def ZtoS(Z,Z0=50):
    Z11=Z[0,0]
    Z12=Z[0,1]
    Z21=Z[1,0]
    Z22=Z[1,1]
    S11 = ( (Z11-Z0)*(Z22+Z0)-Z12*Z21 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S12 = ( 2*Z12*Z0 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S21 = ( 2*Z21*Z0 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S22 = ( (Z11+Z0)*(Z22-Z0)-Z12*Z21 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S=matrix([[S11,S12],[S21,S22]])
    return S


def ABCDtoTransferFct(ABCD,Zs=0,Zl=1e99):
    '''
    To get the transfer function from the ABCD parameters, 
    we can use the equation shown below. In this equation, 
    we consider the impedance from the source side of the network (S) 
    and the load side (L). If the network is terminated to the characteristic 
    impedance on each side, then the two values are equal to the characteristic impedance Z.
    see: https://resources.system-analysis.cadence.com/blog/2020-how-to-calculate-a-transfer-function-from-s-parameters
    
    Parameters
    ----------
    
    ABCD : array 2x2xsize
        phase constant at operating freq 
    Zs : complex
        source impedance
    Zl : complex
        load impedance      
    
    Returns
    -------
    
    complex 
        voltage transfer function H(f) = Vl/Vs
    
    '''
    if shape(ABCD) == (2,2):
        ABCD = array([ABCD])
    A=ABCD[:,0,0]
    B=ABCD[:,0,1]
    C=ABCD[:,1,0]
    D=ABCD[:,1,1]    
    H = Zl / (A*Zl + B + C*Zs*Zl + D*Zs)
    return squeeze(H)

################################################################################
################################################################################
##    AMP DESIGN 
################################################################################
################################################################################

def latexMatrix(a,rnd=None):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :rnd:  rounding digitrs, int
    :returns: LaTeX bmatrix as a string
    """
    set_printoptions(suppress=True)
    if rnd is not None:
        a = around(a,rnd)
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').replace('j','j,').replace('+0.j','').replace('. ','').replace(' 0 ','').replace(' ','').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.rstrip(',').split(',')) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']

    rv =  '\n'.join(rv)
    rv = rv.replace(' ','')
    return rv

### Return a complex type from a number given in magitude and phase (Degrees)
def magphase(A,phi):
    '''Returns a complex number from  magnitude and phase (in degrees)
    '''
    return A*exp(1j*phi*pi/180.0)

### Return a string formatted from a complex in the form Magn /__ Phase deg

################################################################################
def magphase_str(c):
    ''' Returns a nicely formatted string to print complex numbers in ampl. and phase
    '''
    return u'{0:6.3f}\u2220{1:5.1f}\u00B0'.format(abs(c),angle(c)*180/pi)
    
################################################################################    
def magphase_latex(c):
   return r'{0:6.3f}\angle {1:5.1f}^\circ'.format(abs(c),angle(c)*180/pi)     

################################################################################
def magphase_tuple(c):
    ''' Returns a tuple with (magn,phase) to print complex numbers in ampl. and phase
    '''
    return ( abs(c) , angle(c)*180/pi )
    
    
def polar(mag,ang,isDegrees=True):
    '''
    takes a complex number in polar and returns the complex number

    '''
    fac = 1
    if isDegrees:
        fac = pi/180
    return mag*exp(1j*ang*fac)

################################################################################
def splitmatrixarray(S):
    '''
....splits list of matrices into  lists of the individual elements
    currently two by two matirces only
    '''
    S11 = S[:,0,0]
    S12 = S[:,0,1]
    S21 = S[:,1,0]
    S22 = S[:,1,1]
    return S11,S12,S21,S22


### Save touchstone formatted  S-parameter files ###############################
def save_touchstone(filename, flist, slist, annotations = "Touchstone file created by python mwave module "):
    
    '''
    saves a touchstone file of two lists 
       
    Parameters
    ----------
    filename : string 
        name of touchstone file  shall end with .s2p
    flist : array or list
        list of frequency values
    slist : array of 2x2 arrays
        list or array of 2x2 matrices with S-parameters
    annotations : string
        annotations in the header of the file e.g. time 
    Returns 
    -------
    nothing
    
    Examples
    --------
    >>> filename = "touchstone.s2p"
    >>> flist = array([1e9,2e9,2.2e9])
    >>> S1 = array([[0.2,0.3],[0.4,0.5-1j]])
    >>> S2 = array([[0.5,0.333],[0.34,0.35-0.44j]])
    >>> S3 = array([[0.11,0.234],[0.554,0.55-.55j]])
    >>> slist = array([S1,S2,S3])
    >>> save_touchstone(filename, flist, slist)
        
    '''
    
    # check for consistent data 
    if len(flist) != len(slist):
        raise ValueError('length of flist and slist do not match in save touchstone!')
        return
    if shape(slist)[1:3] != (2,2) and ndim(slist) != 1:
        raise ValueError('No 2x2 matrices in touchstone swrite!')
        return
    
    f=open(filename,'w', encoding = "ISO-8859-1")
    noise=False
    f.write('! \n')
    f.write('! Export of Touchstone Data from mwave.py Module Author: S.F. Peik \n')
    f.write('! \n')
    f.write('!'+annotations +'\n')
    f.write('!---------------------------------------------------------------------\n')
    f.write('! symbol freq-unit parameter-type data-format keyword impedance-ohm\n')
    f.write('#        HZ        S              RI          R       50\n')
    f.write('!---------------------------------------------------------------------\n')
    if ndim(slist) == 1:
        # -- One Port parameter -----------
        f.write('! freq       reS11      imS11 \n')
        s11 = slist
        for i in range(len(flist)):
            l = "{:10.1f} {: 3.9f} {: 3.9f}".format(flist[i],real(s11[i]), imag(s11[i]))
            f.write(l+"\n")
    else:
        #--- two-port paramter -----------
        s11,s12,s21,s22 = splitmatrixarray(slist)
        for i in range(len(flist)):
            l = "{:10.1f} {: 3.9f} {: 3.9f} {: 3.9f} {: 3.9f} {: 3.9f} {: 3.9f} {: 3.9f} {: 3.9f} ".format(flist[i],real(s11[i]), imag(s11[i]), real(s12[i]), imag(s12[i]), real(s21[i]), imag(s21[i]), real(s22[i]), imag(s22[i]), )
            f.write(l+"\n")
    f.close()
    return

### Load touchstone formatted  S-parameter files ###############################
def load_touchstone(filename, annotations=False):
    '''
    Loads a touchstone file in two lists 
    
        :filename: Touchstone filename including path (type string)
        
    Returns 
    -------
        tuple with: frequency list (type flaot)  S-matrix list (2x2 Matrix list of S Parameters)
       
    Note
    -----
    
    currently works with 2x2 matrices only
    '''
    
    print("Load Touchstone file ",filename)
    f=open(filename,'r', encoding = "ISO-8859-1")
    noise=False
    #if filename[-2] == '1': 
    #    Twoport = False
    #elif filename[-2] == '2': 
    #    Twoport = True
    #elif filename[-2:] == 'ts': 
    #    Twoport = True
    #else:
    #    print('Load Touchstone: Neither extension s1p or s2p , Exit')
    #    raise NameError('Neither extension s1p or s2p')
    
    try:
        n_ports = int(filename[-2])
    except:
        raise NameError('Not extension sxp or sxp')
    anno = []
    Slist=[];flist=[]
    rad=pi/180.0
    print("Loading ",n_ports,"-Port")
    with open(filename) as fi:
      i = 0
      f = []
      S = []
      line = "!!!!!!!!"
      while len(line)>0:
        line = fi.readline()
        if len(line)<3: continue
        #print(line.strip())
        if line[0]=='!': 
            anno.append(line)
            if line.find('Fmin')>0:
                noise=True
                #print("----- Here Noise Data start ------>")
            continue
        factor = 1.0    
        if line[0]=='#':
            #print("Format is ",line)
            if 'HZ' in line.upper(): factor=1e0
            if 'KHZ' in line.upper(): factor=1e3
            if 'MHZ' in line.upper(): factor=1e6
            if 'GHZ' in line.upper(): factor=1e9
            if 'MA' in line.upper():
                sform ='MA'
            elif 'RI' in line.upper(): 
                sform = 'RI'
            elif 'DB' in line.upper(): 
                sform ='DB'
            else:
                print("Data not in MA or RI Format")
                raise RuntimeError("Data not in MA or RI Format")
                return
            continue
        if len(line) <10: continue ## empty line
        if not(noise): ##### Spara Info
            p=line.split()
            while len(p)-1 < 2*n_ports**2:
                #print("need more Lines", len(p)-1,2*n_ports**2)
                line = fi.readline()
                p.extend(line.split())
            p=[float(x) for x in p]
            #print("f=",p[0],"S11=",p[1], ".....")
            flist.append(float(p[0])*factor)
            # Combine Real Imag ###
            if sform == 'RI':
                # Combine Real Imag into complex number ###   
                _s = array([ p[2*i+1] + 1j*p[2*i+2] for i in range(n_ports**2)])
            elif sform =='DB':
                # Combine dB Phase into complex number ### 
                _s = array([ 10**(p[2*i+1]/20) *  exp(1j*pi/180*p[2*i+2]) for i in range(n_ports**2)])
            elif sform =='MA':
                # Combine dB Phase into complex number ### 
                _s = array([ p[2*i+1] *  exp(1j*pi/180*p[2*i+2]) for i in range(n_ports**2)])
            _S = _s.reshape(n_ports,n_ports)
            Slist.append(_S)
            #print S
        if (noise): ##### Noise Info
            pass
    flist = array(flist)
    Slist = array(Slist)
    if annotations:
        return flist,squeeze(Slist),anno
    else:
        return flist,squeeze(Slist)


#
# mdif load
#
import re
rad = pi/180.0

######################################################################
def mdifbiaslist(filename):
    '''
    Shows the possible bias points of a mdif file
    :param filename: mdif file
    :return: a list of biases
    '''
    f=open(filename,'r')
    line = f.readlines()
    i=0
    biaslist = []
    while i< len(line):
        if 'VAR Vc' in line[i]:
            if not 'Ic' in line[i+1]: 
                raise valueerror('No Vc,Ic VAR defined in mdif')
            valueV = re.findall(r"\d+\.\d+", line[i])[0]
            valueI = line[i+1].rstrip().rstrip("mA").lstrip("VAR Ic=")
            biaslist.append((float(valueV),float(valueI)))
            i += 1   
        i += 1
    if biaslist == []: raise valueerror('No Vc,Ic VAR defined in mdif')
    return biaslist
  
##########Load MDIF Spara #############################################################
def mdifsparlist(filename,Vc,Ic):
    f=open(filename,'r')
    line = f.readlines()
    i=0
    biaslist = []
    while i< len(line):
        if 'VAR Vc' in line[i]:
            try:
                valueV = float(re.findall(r"\d+\.\d+", line[i])[0])
            except:
                valueV = float(re.findall(r"\d+\\d+", line[i])[0])
            try:
                valueI = float(re.findall(r"\d+\.\d+", line[i+1])[0])
            except: 
                valueI = float(re.findall(r"\d+\\d+", line[i+1])[0])
            if float(valueV) == float(Vc) and float(valueI) == float(Ic):
                #print("Biaspoint found", valueV, valueI)
                if not ('BEGIN ACDATA' in line[i+2]): raise ValueError('MDIF Wrong Format no BEGIN ACDATA found ')
                i +=3
                #print(line[i])
                if not '#' in line[i]: raise ValueError('MDIF Wrong Format no # Format found found ')
                if 'HZ'  in line[i]: factor=1e0
                if 'MHZ' in line[i]: factor=1e6
                if 'GHZ' in line[i]: factor=1e9
                if 'MA' in line[i]:
                    sform ='MA'
                elif 'RI' in line[i]: 
                    sform = 'RI'
                else:
                    raise RuntimeError("MDIF Data not in MA or RI Format")
                #print(sform, factor)
                i += 2
                
                ##### Start of spar found reading data ###################
                flist = []
                Slist = []
                while not 'END' in line[i]: 
                    p=line[i].split()
                    p=[float(x) for x in p]
                    #print("f=",p[0],"S11=",p[1], ".....")
                    flist.append(float(p[0])*factor)
                    if sform=='MA':
                        S11=p[1]*exp(1j*p[2]*rad)
                        S21=p[3]*exp(1j*p[4]*rad)
                        S12=p[5]*exp(1j*p[6]*rad)
                        S22=p[7]*exp(1j*p[8]*rad)
                    if sform=='RI':
                        S11=p[1]+p[2]*1j
                        S21=p[3]+p[4]*1j
                        S12=p[5]+p[6]*1j
                        S22=p[7]+p[8]*1j
                    S=matrix([[S11,S12],[S21,S22]])
                    Slist.append(S)
                    i += 1
                return flist, Slist
                ### end of spar data read 
                
        i += 1
    raise ValueError('Specific Vc,Ic not defined in mdif')
    return 

###### Load MDIF Noise #############################################
def mdifnoiselist(filename,Vc,Ic):
    """
    reads and returns the frequency list of noise paramters
    :param filename: filename of mdif file (str)
    :param Vc: Bias Voltage (float)
    :param Ic: Bias Current (float)
    :return: flist, Nfminlist, Gamoptlist, Rnlist
    """
    f=open(filename,'r')
    line = f.readlines()
    i=0
    biaslist = []
    while i< len(line):
        if 'VAR Vc' in line[i]:
            valueV = float(re.findall(r"\d+\.\d+", line[i])[0])
            valueI = line[i + 1].rstrip().rstrip("mA").lstrip("VAR Ic=")

            if float(valueV) == float(Vc) and float(valueI) == float(Ic):
                #print(line[i])
                #print("Biaspoint found", valueV, valueI)
                i+=  2
                while i< len(line):
                    if ('BEGIN NDATA' in line[i]): break
                    i += 1
                if i == len(line): raise ValueError('MDIF no BEGIN NDATA found ')
                i += 1
                if not '#' in line[i]: raise ValueError('MDIF Wrong Format no # Format found found ')
                if 'HZ'  in line[i]: factor=1e0
                if 'MHZ' in line[i]: factor=1e6
                if 'GHZ' in line[i]: factor=1e9
                if 'MA' in line[i]:
                    sform ='MA'
                elif 'RI' in line[i]: 
                    sform = 'RI'
                else:
                    raise RuntimeError("MDIF Data not in MA or RI Format")
            
                i += 2
                ##### Start of spar found reading data ###################
                flist = []
                Nfminlist = []
                Gamoptlist = []
                Rnlist = []
                while not ('END' in line[i]): 
                    p=line[i].split()
                    p=[float(x) for x in p]
                    #print("f=",p[0],"S11=",p[1], ".....")
                    flist.append(p[0]*factor)
                    Nfminlist.append(p[1])    ### min Noisefigure
                    if sform=='MA':
                        Gamoptlist.append(p[2]*exp(1j*p[3]*rad)) ## Gamma Opt. 
                    if sform=='RI':
                        Gamoptlist.append(p[2]+p[3]*1j)
                    Rnlist.append(p[4])
                    i += 1
                if flist == []:
                    print('MDIF: No Noise data defined')
                return flist, Nfminlist, Gamoptlist, Rnlist
                ### end of spar data read 
                
        i += 1
    raise ValueError('Specific Vc,Ic not defined in mdif')
    return 


def ssplit(S):
    '''
    splits list of matrices into  lists of the individual elements
    
    currently two by two matirces only
    
    Examples
    --------
    
    >>> S = array([[[ 0.1+0.j,  0.1+0.j],[ 0.1+0.j, -0.3+0.j]],[[ 0.2+0.j,  0. +0.j],[ 0. +0.j, -0.1+0.j]],[[ 0.3+0.j,  0.4+0.j],[ 0.4+0.j,  0.9+0.j]]]) 
    >>> SSS = ssplit(S)
    >>> for S in SSS: print(S)
    [0.1+0.j 0.2+0.j 0.3+0.j]
    [0.1+0.j 0. +0.j 0.4+0.j]
    [0.1+0.j 0. +0.j 0.4+0.j]
    [-0.3+0.j -0.1+0.j  0.9+0.j]

    '''
    S11 = S[:,0,0]
    S12 = S[:,0,1]
    S21 = S[:,1,0]
    S22 = S[:,1,1]
    return S11,S12,S21,S22

def scombine(s11,s12,s21,s22):
    '''
    combines 4 lists of S-paramters into  1 list with matrices
    
    currently two by two matirces only
    
    Parameters
    ----------
    
    s11 : complex
        S_11 Parameter and so on 
    
    Examples
    --------
    
    >>> S11 = [0.1,0.2,0.3]
    >>> S12 = [0.6,0.7,0.8]
    >>> S21 = [0.1,0.0,0.4]
    >>> S22 = [-0.3,-0.1,0.9]
    >>> scombine(S11,S12,S21,S22)
    array([[[ 0.1+0.j,  0.1+0.j],
            [ 0.1+0.j, -0.3+0.j]],
    <BLANKLINE>
           [[ 0.2+0.j,  0. +0.j],
            [ 0. +0.j, -0.1+0.j]],
    <BLANKLINE>
           [[ 0.3+0.j,  0.4+0.j],
            [ 0.4+0.j,  0.9+0.j]]])

    '''
    S = zeros( (len(s11), 2, 2),dtype="complex" )
    for ii in range(len(s11)):
        S[ii] = array([[s11[ii],s21[ii]],[s21[ii],s22[ii]]])  
    return S

## Smooth data by averaging ########################################################
def smoothTriangle(data, degree):
    triangle=concatenate((arange(degree + 1), arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(sum(point)/sum(triangle))
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed
    
    
### Plot S-Parameter in Cart Plot
def plotspar(flist,Slist=array([0]),funit="MHz",frange=None, phase= False, grid=False, asSmith=False, split = False, smoothing = 1):
    '''
    plot the S-Paramter response quickly 
    asSmith: plot the return losses as smith chart
    smoothing: degree of data smoothing, default 1 no smoothing
    split: split into separate graphs (not implemented)
    '''
    
    flist=array(flist)
    print(Slist.ndim)
    if Slist.ndim == 1:  ### Oneport
        Slist = Slist[:,newaxis,newaxis]
        n_ports = 1
    else:
        n_ports = Slist.shape[1]
    print("# of Ports",n_ports)
    if funit == 'Hz':
        factor = 1.0
    elif funit == 'MHz':
        factor = 1e6
    elif funit == 'GHz':
        factor = 1e9
    elif funit == "1/s":
        factor = 1.0
    else:
        raise ValueError("funit must be one of Hz MHz GHz or 1/s")
    
    ## reduce to selected freq range ####################################
    try:
        if frange is not None:
            idx_min = find_nearest(flist, frange[0])
            idx_max = find_nearest(flist, frange[1])
            flist = flist[idx_min:idx_max] 
            Slist = Slist[idx_min:idx_max]
    except:
        raise ValueError("Invalid Frequency Range selected")

    ## Plot as Smith Charts #############################################
    if asSmith and not grid:
        fig,ax = plt.subplots(figsize=(6,6))
        mysmith = smi.Smith(ax,color='gray')    
        for j in range(n_ports):  # Step thru all S port variations 
            Spar=array([Slist[i][j,j] for i in range(len(Slist))])
            Spar = smoothTriangle(Spar,smoothing)
            ax.plot(real(Spar),imag(Spar),label=f'$S_{{{j+1:}{j+1:}}}$')
        plt.legend(loc=4)
    ## Plot as Grid #####################################################    
    elif grid:
        fig,ax = plt.subplots(n_ports,n_ports,figsize=(10,10))
        for j in range(n_ports):
            for k in range(n_ports):
                ax[j,k].set_title("S"+str(k+1)+str(j+1))
                Spar=array([Slist[i][j,k] for i in range(len(Slist))])
                Spar = smoothTriangle(Spar,smoothing)
                if k != j or (not asSmith):
                    ax[j,k].plot(flist/factor,20*log10(abs(Spar)),color="r")
                    ax[j,k].set_ylim(-15,0)
                    ax[j,k].grid()
                else: 
                    sm = smi.smith(ax[j,k],color="gray")
                    ax[j,k].plot(real(Spar),imag(Spar),color="r")
        plt.tight_layout()
    ## Plot as Cartesian Diagram #######################################    
    else:  # Cartesian Diagram 
        fig,ax = plt.subplots(figsize=(6,6))
        for j in range(n_ports):  # Step thru all S port variations 
            for k in range(n_ports):
                Spar=array([Slist[i][j,k] for i in range(len(Slist))])
                Spar = smoothTriangle(Spar,smoothing)
                if not phase: 
                    ax.plot(flist/factor,20*log10(abs(Spar)),label=f'$S_{{{j+1:}{k+1:}}}$')
                    ax.set_ylabel("dB")
                else:
                    ax.plot(flist/factor,angle(Spar)*180/3.1415,label=f'$\\angle S_{{{j+1:}{k+1:}}}$')
                    ax.set_ylabel("Phase in °")
        plt.legend(loc=4)
        plt.grid()
        ax.set_xlabel("Freq. in "+funit)
    plt.tight_layout()
    return fig,ax
    

def mufactor(S):
    '''
    calculate Mu Factor for an S-parameter Array list
    '''
    
    S11=S[:,0,0];S12=S[:,0,1];S21=S[:,1,0];S22=S[:,1,1]
    Delta=S11*S22-S12*S21
    mu1=(1-abs(S11)**2)/(abs(S22-Delta*conj(S11))+abs(S12*S21))
    mu2=(1-abs(S22)**2)/(abs(S11-Delta*conj(S22))+abs(S12*S21))
    return mu1, mu2

###### determine Stability Circles and mu factor
def AmpStabilityCircle(S,plotit=False): 
    '''
    calculate Stability Circles for given S-Matrix
    return: return Cs,Rs,Cl,Rl,mu1, mu2, fig, ax
            when plotit is false: 
            return Cs,Rs,Cl,Rl,mu1, mu2
    '''
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    #print(u"\n \u25AD\u25AD\u25AD\u25AD Stability \u25AD\u25AD\u25AD\u25AD")
    Delta=S11*S22-S12*S21
    #print("|Delta|=",abs(Delta))
    Cl=conj(S22-Delta*conj(S11))/(abs(S22)**2-abs(Delta)**2)
    Rl=abs(S12*S21/(abs(S22)**2-abs(Delta)**2))
    Cs=conj(S11-Delta*conj(S22))/(abs(S11)**2-abs(Delta)**2)
    Rs=abs(S12*S21/(abs(S11)**2-abs(Delta)**2))
    mu1=(1-abs(S11)**2)/(abs(S22-Delta*conj(S11))+abs(S12*S21));
    mu2=(1-abs(S22)**2)/(abs(S11-Delta*conj(S22))+abs(S12*S21));
    k=(1-abs(S11)**2-abs(S22)**2+Delta**2)/(2*abs(S12*S21));

    fig, ax = (0,0)
    if plotit:
        fig, ax = plt.subplots() 
        plt.tight_layout()
        ax.set_title('Stability Circles') 
        fig.set_facecolor('white')
        Z0=1 
        mysmith=smi.smith(ax,'smith',Z0,0.5)
        mysmith.addcircle(Cl,Rl)
        mysmith.addcircle(Cs,Rs,'r')
        #plt.savefig('stabcircles.pdf')
        return Cs,Rs,Cl,Rl,mu1, mu2, fig, ax
    return Cs,Rs,Cl,Rl,mu1, mu2


##### determine Noise Circle and Noise Number N
def AmpNoiseCircle(S,FmindB,Gamopt,rn,FsetdB,plotit): 
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    Delta=S11*S22-S12*S21
    F=10**(FsetdB/10)
    Fmin=10**(FmindB/10)
    N=(F-Fmin)/4/rn*abs(1+Gamopt)**2
    Cf=Gamopt/(N+1)
    Rf=sqrt(N*(N+1-abs(Gamopt)**2))/(N+1)
    fig, ax = (0,0)
    if plotit:
        fig, ax = plt.subplots() 
        plt.tight_layout()
        ax.set_title('Noise Circle', fontsize=15) 
        fig.set_facecolor('white')
        Z0=1 
        mysmith=smi.smith(ax,'smith',Z0,0.5)
        mysmith.addcircle(Cf,Rf)
    return Cf, Rf, N, fig, ax

##### determine Gain Circle 
def AmpGainCircleSource(S,GaindB): 
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    if S12 != 0:
        raise ValueError('Device is not Unilateral!')
        return
    Delta=S11*S22-S12*S21
    Gain=10**(GaindB/10);
    Gsmax=1/(1-abs(S11)**2);
    GsmaxdB=10*log10(Gsmax);
    gs=10**(GaindB/10)/Gsmax;
    if 1-gs <0:
        return nan,nan,nan
    Cs=gs*conj(S11)/(1-(1-gs)*abs(S11)**2);
    Rs=sqrt(1-gs)*(1-abs(S11**2))/(1-(1-gs)*abs(S11)**2);
    return Cs, Rs, GsmaxdB

##### determine Gain Circle 
def AmpGainCircleLoad(S,GaindB): 
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    if S12 != 0:
        raise ValueError('Device is not Unilateral!')
        return
    Delta=S11*S22-S12*S21
    Gain=10**(GaindB/10);
    Gsmax=1/(1-abs(S22)**2);
    GsmaxdB=10*log10(Gsmax);
    gs=10**(GaindB/10)/Gsmax;
    if 1-gs <0:
        return nan,nan, nan
    Cs=gs*conj(S22)/(1-(1-gs)*abs(S22)**2);
    Rs=sqrt(1-gs)*(1-abs(S22**2))/(1-(1-gs)*abs(S22)**2);
    return Cs, Rs, GsmaxdB


def transducerGain(S, Gams, Gaml):
    '''
    Calculates the transducer gain for a given S-Matrix
    and Γ_s and Γ_l

    Parameters
    ----------
    :param S: (2x2 Matrix complex) or array of matrices
    :param Gams: Γ_s (complex)
    :param Gaml: Γ_l (complex)
    :return: Transducer Gain Gt in dB (float or array float)

    Example
    -------
    >>> S = matrix([[0.3+0.2j,0.02],[4.0+2.6j,0.7j]])
    >>> Gam_s, Gam_l = 0.3, 0.1-0.5j
    >>> GtdB = transducerGain(S, Gam_s, Gam_l)
    >>> GtdB
    16.51983047097562

    '''
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    Gamin = S11 + (S12 * S21 * Gaml) / (1 - S22 * Gaml)
    Gamout = S22 + (S12 * S21 * Gams) / (1 - S11 * Gams)
    Gs = (1 - abs(Gams) ** 2) / abs(1 - Gamin * Gams) ** 2
    G0 = abs(S21) ** 2
    Gl = (1 - abs(Gaml) ** 2) / abs(1 - S22 * Gaml) ** 2
    Gt = Gs * G0 * Gl
    GsdB = 10 * log10(Gs)
    G0dB = 10 * log10(G0)
    GldB = 10 * log10(Gl)
    GtdB = GsdB + G0dB + GldB
    return GtdB

####### Bilateral Max. Gain Design
def AmpMaxgain(S, verbose = False):
    r'''
    Calculates Maximum Gain input and output loads
    
    Parameters
    ----------
    
    S : 2x2 Matrix
        S-Matrix containing linear transistor S data 
    
    Returns
    -------
    tuple 
        :math:`\Gamma_s` (type complex), :math:`\Gamma_l` (type complex), GtdB (type float)
    
    Returns (0,0,0) if no stable solution is found
    Examples
    --------
    
    >>> S = matrix([[0.3+0.2j,0.02],[4.0+2.6j,0.7j]])
    >>> Gams,Gaml,Gmax = AmpMaxgain(S)
    >>> print(around(Gmax,2))
    17.21
    >>> print(magphase_str(Gams),magphase_str(Gl)) 
     0.399∠-13.7°  0.718∠-86.8°
    '''
 
    S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
    Delta=S11*S22-S12*S21
    if verbose: print(u"\n\u25AD\u25AD\u25AD\u25AD Conj. Matching \u25AD\u25AD\u25AD\u25AD\n")
    B1=1+abs(S11)**2-abs(S22)**2-abs(Delta)**2
    B2=1+abs(S22)**2-abs(S11)**2-abs(Delta)**2
    C1=S11-Delta*conj(S22)
    C2=S22-Delta*conj(S11)
    stablegain=0
    if verbose: print(B1,B2,C1,C2,Delta)
    if verbose: print("Solution 1:")
    Gams1=(B1-sqrt(B1**2-4*abs(C1)**2+0j))/2/C1
    Gaml1=(B2-sqrt(B2**2-4*abs(C2)**2+0j))/2/C2
    if verbose: print(u"   \u0393s1=",magphase_str(Gams1))
    if verbose: print(u"   \u0393l1=",magphase_str(Gaml1))
    if verbose: print("Solution 2:")
    Gams2=(B1+sqrt(B1**2-4*abs(C1)**2+0j))/2/C1
    Gaml2=(B2+sqrt(B2**2-4*abs(C2)**2+0j))/2/C2
    if verbose: print(u"   \u0393s2=",magphase_str(Gams2))
    if verbose: print(u"   \u0393l2=",magphase_str(Gaml2))
    if (abs(Gams1)<0.99 and abs(Gaml1)<0.99):
        if verbose: print(">> Choosen Sol.1:")
        if verbose: print( u"   \u0393s=",magphase_str(Gams1))
        if verbose: print( u"   \u0393l=",magphase_str(Gaml1))
        Gaml=Gaml1
        Gams=Gams1
        stablegain=1;
    if (abs(Gams2)<0.99 and abs(Gaml2)<0.99):
        if verbose: print(">> Choosen Sol.2:")
        if verbose: print(u"   \u0393s=",magphase_str(Gams2))
        if verbose: print(u"   \u0393l=",magphase_str(Gaml2))
        Gaml=Gaml2
        Gams=Gams2
        stablegain=1
    if not(stablegain):
        if verbose: print(" ***** Geht Nicht Unstable *****")
        raise ValueError("Unstable")
        return 0,0,0
    else:
        if verbose: print("\n=== Transducer Power Gain: ===")
        Gamin=S11+(S12*S21*Gaml)/(1-S22*Gaml)
        Gamout=S22+(S12*S21*Gams)/(1-S11*Gams)
        Gs=(1-abs(Gams)**2)/abs(1-Gamin*Gams)**2
        G0=abs(S21)**2
        Gl=(1-abs(Gaml)**2)/abs(1-S22*Gaml)**2
        Gt=Gs*G0*Gl
        GsdB=10*log10(Gs)
        G0dB=10*log10(G0)
        GldB=10*log10(Gl)
        GtdB=GsdB+G0dB+GldB
        if verbose: print(u"Gt={0:.2f}*{1:.2f}*{2:.2f}={3:.2f}".format(Gs,G0,Gl,Gt))
        if verbose: print(u"GtdB={0:.2f}dB+{1:.2f}dB+{2:.2f}dB={3:.2f}dB".format(GsdB,G0dB,GldB,GtdB) )
        return Gams,Gaml,GtdB
    

#### Matching Circuits ##############

def LMatching(Zl,Z0=50, equaltype=False):
    '''
    see: Whites, EE 481/581, Lecture 7: Transmission Line Matching Using Lumped L Networks
    returns a dictionary with possible solutions
    when equaltype is True, both elements are chosen of same type, e.g. capacitor and capacitor
    '''
    Rl = real(Zl)
    Xl = imag(Zl)
    Yl = 1/Zl
    Gl = real(Yl)
    Bl = imag(Yl)
    if Rl > Z0 or equaltype:
        B1 = ( Xl + sqrt(Rl/Z0) * sqrt(Rl*(Rl-Z0)+Xl**2) ) / (Rl**2 + Xl**2)
        X1 = 1/B1 + (Xl*Z0)/Rl - Z0/(B1*Rl)
        B2 = ( Xl - sqrt(Rl/Z0) * sqrt(Rl*(Rl-Z0)+Xl**2) ) / (Rl**2 + Xl**2)
        X2 = 1/B2 + (Xl*Z0)/Rl - Z0/(B2*Rl)  
        if isnan(B1): raise ValueError("No Solution")
        type = 'shunt-series' 
    else:
        X1 = + sqrt(Rl*(Z0-Rl)) - Xl
        B1 = + 1/Z0 * sqrt((Z0-Rl)/Rl)
        X2 = - sqrt(Rl*(Z0-Rl)) - Xl
        B2 = - 1/Z0 * sqrt((Z0-Rl)/Rl)
        if isnan(B1): ValueError("No Solution")
        type = 'series-shunt'
    return {'sol1':(X1,B1),'sol2':(X2,B2)}, type


def AmpStubmatching(Gammamatch,plotit=False):
    r'''
    Performs a complete open stub - line matching for a given desired input :math:`\Gamma`
    
    Plots the smith chart with the constructed matching network 
    
    Parameters
    ----------
    
    Gammamatch : complex
        desired input :math:`\Gamma` 
    plotit : Boolean
        Flag for plotting, if True a Plot is generated
    
    Returns
    -------
    
    tuple 
        line length of stub (float), line length of line (float), fig, ax (handles of smith plot)
    '''

    Z0=1
    Z1=Z0
    ### finding length of open stub line
    for len1 in arange(0.001,1,0.001):
        betal=2*pi*len1
        Zstub=-1j*Z0/tan(betal)
        Z2=Z1*Zstub/(Z1+Zstub)
        Gam=(Z2-Z0)/(Z2+Z0)
        if abs(Gam)>abs(Gammamatch): break
    ### finding length of inserted line
    if(angle(Gam)>angle(Gammamatch)):
        betal=(angle(Gam)-angle(Gammamatch))/2
    else:
        betal=(angle(Gam)-angle(Gammamatch))/2+pi
    len2=round(betal/2/pi,3)
    #print(u"stub length={0:.3f} \u03BB inserted line length={1:.3f}\u03BB".format(len1,len2))
    if plotit:
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)
        plt.tight_layout()
        fig.set_facecolor('white')
        plt.clf()
        ax = plt.axes() 
        Z0=1 
        mysmith=smi.smith(ax,'both',Z0) 
        #mysmith.addanglering()
        Zl=mysmith.addstart(1) 
        Z2=mysmith.addstubopen(1,len1,1)    
        Z3=mysmith.addline(Z2,len2,1)
        mysmith.addpoint(1,'$Z_{50}$','NE')
        mysmith.addpoint(Z2,'$Z_2$','NE') 
        mysmith.addpoint(Z3,'$Z_3$','NE')
        return len1,len2, fig, ax
    else:
        return len1,len2,0,0


############################################################################################
## Common-Differential Tools
############################################################################################


def impedanceFromS_series(S,Z0=50):
    '''
    Find series impedance from 2-Port S-Paramter
    See Cispr 17 6.3.1.3 for more information
    
    S: 2x2  Matrix array S-Matrix
    Z0: Reference impedance (float)
    returns: series impedance in Ohms
    '''
    R = (S[:,0,0] + S[:,1,1]) / 2
    T = (S[:,0,1] + S[:,1,0]) / 2
    Z = Z0 * ((1+R)**2 - T**2) / (2*T) # Reihenschaltung
    return Z


def impedanceFromS_shunt(S,Z0=50):
    '''
    Find shunt impedance from 2-Port S-Paramter
    See Cispr 17 6.3.1.3 for more information
    
    S: 2x2  Matrix array S-Matrix
    Z0: Reference impedance (float)
    returns: shunt impedance in Ohms
    '''
    R = (S[:,0,0] + S[:,1,1]) / 2
    T = (S[:,0,1] + S[:,1,0]) / 2
    Zx = Z0 * (2*T) / ((1-R)**2-T**2)  # Shuntschaltung 
    return Z


############################################################################################


def fourPortS_Matrix_from_TwoPortMatrix(S12_file,S13_file,S14_file):
    '''
    assemble a 4-Port S-Matrix from three seperate 2-Port S-Matrices from a **symmetric** network 
     1 o----|     |---o 2
            |  S  |
     3 o----|     |---o 4
     
     parameter: 
     S12_file: filename of first S-parameter Touchstone file measured Port 1 and Port 2
     S13_file: filename of first S-parameter Touchstone file measured Port 1 and Port 3
     S14_file: filename of first S-parameter Touchstone file measured Port 1 and Port 4
     
     returns: 
     array (1dim): frequency vector 
     array (4dim): 4-Port S-Matrix array over frequency
     
    '''
    fm,Ss12 = load_touchstone(S12_file)
    fm,Ss13 = load_touchstone(S13_file)
    fm,Ss14 = load_touchstone(S14_file)

    S11 = Ss12[:,0,0]
    S12 = Ss12[:,0,1]
    S21 = Ss12[:,1,0]
    S22 = Ss12[:,1,1]

    S11 = Ss13[:,0,0]
    S13 = Ss13[:,0,1]
    S31 = Ss13[:,1,0]
    S33 = Ss13[:,1,1]

    S22 = Ss13[:,0,0]
    S24 = Ss13[:,0,1]
    S42 = Ss13[:,1,0]
    S44 = Ss13[:,1,1]

    S33 = Ss12[:,0,0]
    S34 = Ss12[:,0,1]
    S43 = Ss12[:,1,0]
    S44 = Ss12[:,1,1]

    S22 = Ss14[:,0,0]
    S23 = Ss14[:,0,1]
    S32 = Ss14[:,1,0]
    S33 = Ss14[:,1,1]

    S11 = Ss14[:,0,0]
    S14 = Ss14[:,0,1]
    S41 = Ss14[:,1,0]
    S44 = Ss14[:,1,1]


    S =    array([[S11,S12,S13,S14],
                  [S21,S22,S23,S24],
                  [S31,S32,S33,S34],
                  [S41,S42,S43,S44]]).transpose(2,0,1)
    return fm,S
                  


def threePortS_Matrix_from_TwoPortMatrix(S12_file,S13_file,S23_file):
    '''
    assemble a 4-Port S-Matrix from three seperate 2-Port S-Matrices from a symmetric network 
            |     |---o 2
     1 o----|  S  |
            |     |---o 3
     
     parameter: 
     S12_file: filename of first S-parameter Touchstone file measured Port 1 and Port 2
     S13_file: filename of first S-parameter Touchstone file measured Port 1 and Port 3
     S23_file: filename of first S-parameter Touchstone file measured Port 2 and Port 3
     returns: 
     array (1dim): frequency vector 
     array (4dim): 3-Port S-Matrix array over frequency
     
    '''
    fm,Ss12 = load_touchstone(S12_file)
    fm,Ss13 = load_touchstone(S13_file)
    fm,Ss23 = load_touchstone(S23_file)

    S11 = Ss12[:,0,0]
    S12 = Ss12[:,0,1]
    S21 = Ss12[:,1,0]
    S22 = Ss12[:,1,1]

    #S11 = Ss13[:,0,0]
    S13 = Ss13[:,0,1]
    S31 = Ss13[:,1,0]
    S33 = Ss13[:,1,1]

    #S22 = Ss23[:,0,0]
    S23 = Ss23[:,0,1]
    S32 = Ss23[:,1,0]
    #S33 = Ss23[:,1,1]
    # In comments redundant measurements
   

    S =    array([[S11,S12,S13],
                  [S21,S22,S23],
                  [S31,S32,S33]]).transpose(2,0,1)
    return fm,S
                


############################################################################################
def singleEndedToMixedMode(S):
    '''
    converts a single ended S-Parameter File to mixed Mode S-Parameter
    S: single Ended 4-Port S-Parameter File
    retruns: 
    Smixed s. lit. 

    1 o----|     |---o 2
              S 
    3 o----|     |---o 4
    
    returns: 
    Scc: common-common S matrix
    Sdd: diff-diff S matrix
    Scd: common-diff S matrix
    Sdc: diff-common S matrix
    
    see also https://coppermountaintech.com/wp-content/uploads/2022/06/BalancedMeas.pdf
    '''
    
    if S.ndim ==2:
        S = array([S])
    if S.ndim != 3:
        raise ValueError("Matrix must be n_freq x 4 x 4 here " + str(shape(S)))
    if not shape(S)[1:3] == (4,4):
        raise ValueError("Matrix must be n_freq x 4 x 4 here "+str(shape(S)))
    
    S11,S12,S13,S14 = S[:,0,0], S[:,0,1], S[:,0,2], S[:,0,3]
    S21,S22,S23,S24 = S[:,1,0], S[:,1,1], S[:,1,2], S[:,1,3]
    S31,S32,S33,S34 = S[:,2,0], S[:,2,1], S[:,2,2], S[:,2,3]
    S41,S42,S43,S44 = S[:,3,0], S[:,3,1], S[:,3,2], S[:,3,3]
    
    SDD11=0.5*(S11-S13-S31+S33)
    SDD22=0.5*(S22-S24-S42+S44)
    SDD21=0.5*(S21-S23-S41+S43)
    SDD12=0.5*(S12-S14-S32+S34)
    SCD11=0.5*(S11-S13+S31-S33)
    SCD22=0.5*(S22-S24+S42-S44)
    SCD21=0.5*(S21-S23+S41-S43)
    SCD12=0.5*(S12-S14+S32-S34)
    
    SDC11=0.5*(S11+S13-S31-S33)
    SDC22=0.5*(S22+S24-S42-S44)
    SDC21=0.5*(S21+S23-S41-S43)
    SDC12=0.5*(S12+S14-S32-S34)
    SCC11=0.5*(S11+S13+S31+S33)
    SCC22=0.5*(S22+S24+S42+S44)
    SCC21=0.5*(S21+S23+S41+S43)
    SCC12=0.5*(S12+S14+S32+S34)
    
    Scc = [ [SCC11,SCC12],[ SCC21,SCC22] ] 
    Sdd = [ [SDD11,SDD12],[ SDD21,SDD22] ]  
    Scd = [ [SCD11,SCD12],[ SCD21,SCD22] ] 
    Sdc = [ [SDC11,SDC12],[ SDC21,SDC22] ] 
    Scc = transpose(Scc,(2,0,1))
    Sdd = transpose(Sdd,(2,0,1))
    Sdc = transpose(Sdc,(2,0,1))
    Scd = transpose(Scd,(2,0,1))

    return squeeze(Scc),squeeze(Sdd),squeeze(Scd),squeeze(Sdc)


#####################################################################################################


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
 
