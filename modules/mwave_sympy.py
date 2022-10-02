"""

.. module:: Useful Microwave Tools Module 
   :platform: Unix, Windows
   :synopsis: A toolbox of MW-Tools
.. moduleauthor:: Soeren Peik<speik@hs-bremen.de>

"""

from numpy import  sqrt,pi,log,matrix
from sympy import Matrix, Rational
import types

#### some constants

mu0  = 4*pi*1e-7    # Permeability of free space
c    = 299792458.0  # velocity of light
eps0 = 1/mu0/c/c    # Permittivity of free space
k    = 1.3806488e-23 # Boltzman Constant
eta0 = sqrt(mu0/eps0) # Free Space impedance
 

#
# calculates microstrip line impedance,e_eff from Wheeler formula
# S. Peik
#
# Oct 2013
#


###########################################################################
def msimpedance(w,h,er):
  '''
  calculates microstrip line impedance,e_eff from Wheeler formula 

    :param w: width of MS-line
    :type w: float
    :param h: height of substrate
    :type h: float
    :param er: Epsilon relative of substrate
    :type er: float
    :returns:  Z0, eps_eff 
    :raises: non

  '''

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


### ABCD Matrix for Series Element #########################
def ABCDseries(Z):
    return matrix([[1,Z],[0,1]])

### ABCD Matrix for Shunt Element #########################
def ABCDshunt(Y):
    return matrix([[1,0],[Y,1]])

############################################################
def cascade(ABCDlist):
    p = matrix([[1,0],[0,1]])
    for ABCD in ABCDlist:
        p = p * ABCD
    return p

############################################################
def ABCDJInverter(J):
    return(matrix([[0,1/(1j*J)],[-1j*J,0]]))

def ABCDKInverter(K):
    return(matrix([[0,1j*K],[-1/1j/K,0]]))

############################################################
def ABCDtoS(ABCD,Z0):
    if isinstance(ABCD[1,0], types.FunctionType):
        print("Function included")
        return
    else:
	    pass
    A=ABCD[0,0]
    B=ABCD[0,1]
    C=ABCD[1,0]
    D=ABCD[1,1]
    S11=(A+B/Z0-C*Z0-D) /(A+B/Z0+C*Z0+D)
    S12=2*(A*D-B*C)     /(A+B/Z0+C*Z0+D)
    S21=2               /(A+B/Z0+C*Z0+D)
    S22=(-A+B/Z0-C*Z0+D)/(A+B/Z0+C*Z0+D)
    S=Matrix([[S11,S12],[S21,S22]])
    return S

############################################################
def StoABCD(S,Z0):
    S11=S[0,0]
    S12=S[0,1]
    S21=S[1,0]
    S22=S[1,1]
    A =      ((1+S11)*(1-S22)+S12*S21) / 2/S21
    B = Z0*  ((1+S11)*(1+S22)-S12*S21) / 2/S21
    C =      ((1-S11)*(1-S22)-S12*S21) / 2/S21  / Z0
    D =      ((1-S11)*(1+S22)+S12*S21) / 2/S21
    ABCD=Matrix([[A,B],[C,D]])
    return ABCD



#############################################################
def ZtoS(Z,Z0):
    Z11=Z[0,0]
    Z12=Z[0,1]
    Z21=Z[1,0]
    Z22=Z[1,1]
    S11 = ( (Z11-Z0)*(Z22+Z0)-Z12*Z21 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S12 = ( 2*Z12*Z0 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S21 = ( 2*Z21*Z0 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S22 = ( (Z11+Z0)*(Z22-Z0)-Z12*Z21 ) / ((Z11+Z0)*(Z22+Z0)-Z12*Z21)
    S=Matrix([[S11,S12],[S21,S22]])
    return S





