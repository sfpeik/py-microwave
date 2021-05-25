# -*- coding: utf-8 -*-
'''
S. Peik Aug 2017
V1.2 Dez 2017: include M-Extract
'''
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.optimize import least_squares
from pprint import pprint

import schemdraw as schem # Schemdraw packagw  install via: sudo pip3 install schemdraw 
import schemdraw.elements as e 


###############################################################################
def coth(x):
    ''' calculate the coth '''
    return 1/np.tanh(x)


###############################################################################
def gcoeffs(n, LdB):
    r'''Calculates the g-coefficients for doubly terminated filters, 
    Chebychev or Butterworth is supported
    
    Parameters
    ----------
    n : integer
        Order of Filter
    LdB : float
        Ripple of Filter in dB, 0 means Butterworth

    Returns 
    -------
    array of floats
        g-coefficients
    
    Examples
    --------
        
    >>> g = gcoeffs(5,0.5)
    >>> print(np.around(g,4))
    [1.     1.7058 1.2296 2.5409 1.2296 1.7058 1.    ]
    >>> for ii in range(1,8):
    ...     g = gcoeffs(ii,3.0)
    ...     print("Order ",ii," :", np.around(g,4))
    Order  1  : [1.     1.9954 1.    ]
    Order  2  : [1.     3.1014 0.5339 5.8095]
    Order  3  : [1.     3.3489 0.7117 3.3489 1.    ]
    Order  4  : [1.     3.4391 0.7483 4.3473 0.592  5.8095]
    Order  5  : [1.     3.4815 0.7619 4.5378 0.7619 3.4815 1.    ]
    Order  6  : [1.     3.5047 0.7685 4.6063 0.7929 4.4643 0.6033 5.8095]
    Order  7  : [1.     3.5187 0.7722 4.6392 0.8038 4.6392 0.7722 3.5187 1.    ]

    '''


    g = np.zeros(n + 2)
    g[0] = 1.
    if (LdB == 0):   # Butterworth
        g = [2. * np.sin(((2. * k - 1) * np.pi) / (2. * n)) for k in range(0, n + 1)]
        g[0] = 1.
        g.append(1.0)
    else:  # Chebychev
        beta = np.log(coth(LdB / 17.37))
        gamma = np.sinh(beta / 2 / n)
        a = [np.sin(((2. * k - 1) * np.pi) / (2. * n)) for k in range(0, n + 1)]
        b = [gamma ** 2 + (np.sin(k * np.pi / n)) ** 2 for k in range(0, n + 1)]
        g[1] = 2 * a[1] / gamma
        for k in range(2, n + 1):
            g[k] = 4 * a[k - 1] * a[k] / b[k - 1] / g[k - 1]
        g[n + 1] = 1
        if (n % 2 == 0):
            g[n + 1] = (coth(beta / 4)) ** 2
    return g


###############################################################################
def MofChebychev(n,rippledB):
    '''
    A Quick way to generate a Chebychev all-pole coupling matrix
     
    Parameters
    ----------
    n : integer
        Order of Filter
    rippledB : float
        Ripple of Filter in dB, 0 means Butterworth
        
    Returns 
    -------
    n+2 x n+2 array 
        coupling matrix
    
    Examples
    --------
    >>> M = MofChebychev(3,1.0)
    >>> np.around(M,3)
    array([[0.   , 0.703, 0.   , 0.   , 0.   ],
           [0.703, 0.   , 0.705, 0.   , 0.   ],
           [0.   , 0.705, 0.   , 0.705, 0.   ],
           [0.   , 0.   , 0.705, 0.   , 0.703],
           [0.   , 0.   , 0.   , 0.703, 0.   ]])
    '''  
    
    ### Determine G-Values for given Prototype ################
    g = gcoeffs(n, rippledB)
    ### Setup Coupling Matrix from G-Values ###################
    M = np.matrix(np.zeros((n+2, n+2)),dtype = float)
    for k in range(0, n+1):
        M[k, k+1] = 1 / np.sqrt(g[k] * g[k + 1])
        M[k+1, k] = M[k,k+1]
    return M
    

###############################################################################
def RespM2(M, f0, BW, Q , frange, NRNlist=[]):
    '''Calculating Response of the N+2 Matrix using Cameron, Kudsia, Mansour, Including NRNs
     
    Parameters
    ----------
        M : 2-dim array of floats
            Coupling Matrix n+2 x n+2
        f0 : float
            centre frequency
        BW: float
            Realtive Bandwidth of the filter
        Q:  float
            Quality factor of resonators, currently identical for all resonators
        frange: array of floats 
            vector with  frequency samples to be analysed
        NRNlist: list
            list of nonresonating nodes (excl. source, load), e.g. NRNlist = [2,3,5]
            where 0 is source and -1 (order+1) is load
        
    Returns 
    -------
    array of complex
        S-paramter-matrix array of same length as frange
      
    Examples
    --------
    >>> M = MofChebychev(5,1.0)
    >>> frange = np.arange(1.8,2.3,0.1)
    >>> S = RespM2(M, 2.0, 0.1, 1e9, frange)
    >>> np.around(S,3)
    array([[[-0.884+0.467j, -0.002-0.004j],
            [-0.002-0.004j, -0.884+0.467j]],
    <BLANKLINE>
           [[ 0.368+0.549j, -0.623+0.418j],
            [-0.623+0.418j,  0.368+0.549j]],
    <BLANKLINE>
           [[-0.   +0.j   , -1.   -0.j   ],
            [-1.   -0.j   , -0.   +0.j   ]],
    <BLANKLINE>
           [[ 0.213-0.081j, -0.348-0.91j ],
            [-0.348-0.91j ,  0.213-0.081j]],
    <BLANKLINE>
           [[-0.853-0.522j, -0.004+0.006j],
            [-0.004+0.006j, -0.853-0.522j]]])
    
    >>> NRNs = [1,5]
    >>> S = RespM2(M, 2.0, 0.1, 1e9, frange, NRNs)
    >>> S[4] # 5th freq. point only
    array([[0.56308884+0.82294195j, 0.06229444-0.04262428j],
           [0.06229444-0.04262428j, 0.56308884+0.82294195j]])
    '''

    n, n2 = M.shape
    #print ("Order=",n-2)
    if n!=n2: 
        print("ERROR: M is not Square")
        return(1)
    for row in range(n):
        for col in range(n):
            if M[row,col] != M[col,row]:
               print("ERROR: M is not symmetric about diagonal")
               return(1)
    Rs = 1 
    Rl = 1 
    ######## Q Vector Setup ################

    if isinstance(Q, (list, tuple, np.ndarray)):
         if len(Q) != n:
             print("Q Vector length ",len(Q)," is not equal to order ",n,)
             exit(1)
         Qvec = np.array(Q)
    else:
        Qvec = np.ones(n)*Q
    #print("Q-Vector=",Qvec)
    
    #%%%%%%% Matrix fill %%%%%%%%%%%%%%%%%%%
    I = np.eye(n,dtype=complex)
    I[0,0] = 0
    I[-1,-1] = 0
    ######### consider non resonating nodes 
    for nrn in NRNlist: 
        #print("Insert NRN at Res",nrn)
        I[nrn,nrn] = 0
    ######### Resistor Load Matrix
  
    J = 0*np.identity(n,dtype=complex) 
    J[0,0] = 1j*Rs
    J[-1,-1] = 1j*Rl
    
    ## Add losses to resonators due to Q ############# see http://www.jpier.org/PIER/pier115/18.11021604.pdf
    G = 0*np.identity(n,dtype=complex) 
    for res in range(1,n-1):
        if f0 == 0:
            G[res,res] = 1/BW/Qvec[res]
        else:
            G[res,res] = 1/BW/Qvec[res]
    Mpr = M - 1j*G
    #%%%%%%% Sweep frequencies for S-Paramter %%%%%%%%
    flist = []
    s21list = np.zeros(len(frange),dtype="complex")
    s11list = np.zeros(len(frange),dtype="complex")
    s22list = np.zeros(len(frange),dtype="complex")
    Slist = np.zeros( (len(frange), 2, 2),dtype="complex" )
    ii = 0
    for f in frange:
        if f0 == 0:
            lam = 1 / BW * f
        else:
            lam = 1 / BW * (f / f0 - f0 / f)   
        
        Ainv = np.linalg.inv( lam*I - J + Mpr)
        S11 = 1 + 2j * Rs * Ainv[0,0]       
        S22 = 1 + 2j * Rs * Ainv[-1,-1]    
        S21 = -2j * np.sqrt(Rs * Rl)  * Ainv[-1,0]
        s21list[ii] = S21
        s11list[ii] = S11
        s22list[ii] = S22
        Slist[ii] = np.matrix([[S11,S21],[S21,S22]])  
        ii +=1
        
    return Slist

########################################################################
def Mastable(M,rou=4):
    r'''
    Create a LaTeX table of a coupling matrix as defined above
    
    Parameters
    ----------
    M : 2-dim array of floats
        Coupling Matrix n+2 x n+2
    rou: integer
        round of after rou digits, default  is 4 

    Returns 
    -------
    string
        Latex code of the CM table
      
    Examples
    --------
    >>> M = MofChebychev(4,0.5)
    >>> M[1,4] = -0.17
    >>> M[4,1] = M[1,4]
    >>> frange = np.arange(3,5,0.01)
    >>> S = RespM2(M, 4.0, 0.1,1e9,frange)
    >>> Lcode = Mastable(M)
    >>> plt.rc('text', usetex=True)
    >>> fig,(ax1,ax2) = plt.subplots(1,2, figsize=(12,6))
    >>> tt = ax1.plot(frange,20*np.log10(abs(S[:,0,0])))
    >>> tt = ax1.plot(frange,20*np.log10(abs(S[:,1,0])))
    >>> tt = ax1.grid()
    >>> tt = ax1.set_ylim(-50,2)
    >>> tt = ax2.axis('off')
    >>> tt = ax2.text(-0.1,0.15,Lcode,fontsize=12)
    >>> fig.savefig('latexplot.png',dpi=300)
    
    .. figure:: latexplot.png
    
    '''
    
    n,m = np.shape(M)
    ttt = r"\sf \renewcommand{\arraystretch}{2} " 
    ttt +=r"\begin{tabular}{"+"c|"*(n+1)+"} "
    ttt += r" & \bf S"
    for c in range(1,m-1):
        ttt += r" & \bf "+str(c)
    ttt += r" & \bf  L  \\ \hline "
    for row in range(n):
        if row == 0: ttt += r" \bf S & "
        elif row == n-1: ttt += r" \bf L & "
        else: ttt += r"\bf " + str(row)+" & "
        for col in range(m):
            ttt += str(np.around(M[row,col],rou)) +" &"
        ttt = ttt.strip("&")
        ttt += r" \\ \hline "
    ttt += r" \end{tabular}"
    return ttt
    
    
########################################################################
def fit_to_M(flist,Smeas,Minitial,f0,BW,Qinitial,NRNlist,Mindiceslist,is_symmetric=False,deci = 1,fstart=0,fstop=0,maxerr = 5):
    '''
    Tries to fit a coupling matrix to a given S-parameter response using least-square fit
    
    Parameters
    ----------
    
        Mindices : n+2x2+2 array
            The indices of the to be opt M-parameter M(i,j)
                  always symmetric, i.e. i<=j
                  all other M are untouched
        
        flist : numpy array float
            freq points vector, 
        
        Smeas :  numpy 2x2xsamples array complex
            measured Sparameter ,
        
        Minitial: 
            first Guess of Coupling matrix size n+2, numpy array (n+2,n+2), n is order
        
        f0: centre freq, float
        
        BW: Bandwidth, float
        
        Qinitial: initial guess of Q factor
        
        NRNlist, list of non resonatong nodes, python list
        
        Mindiceslist: list of M-parameters to vary, list of 2-tuples, e.g. ((1,1),(1,2),(2,2))
        
        is_symmetric: set True if the filter is symmetric in the port1-port2 layout
        
        deci: decimation of sampled data for faster processing
        
        fstart: start freq of matching response
        
        fstop:
            stop freq of matching response 

    Returns 
    -------
        complex matrix, float, float, float, object :
            M,Q, phaseoffsetinput, phaseoffsetoutput, res (leastsq optimisation information see scipy.least_square)
    Examples
    --------
    
    >>> ## to verify the function we generate a simulated measurement of a 4 pole filter including noise
    >>> frange = np.arange(1.6e9,2.4e9,0.002e9)
    >>> ## assume some invented coupling matrix ############
    >>> n = 3 # filter order
    >>> Mmeas = np.array( [[0.   , 0.6   , 0.   , 0.    , 0.   ], \
                           [0.6  , -0.14 , 0.705, 0.    , 0.   ], \
                           [0.   , 0.705 , 0.06 , 0.505 , 0.   ], \
                           [0.   , 0.    , 0.505, -0.211, 0.703], \
                           [0.   , 0.    , 0.   , 0.703 , 0.   ]])  
    >>> ### Calculate and plot response ####################
    >>> Smeas = RespM2(Mmeas, 2.0e9, 0.1, 1e9, frange)
    >>> ###add phase offsets due to lines ##################
    >>> line1_l = 0.03 # line 1 length in m 
    >>> line2_l = 0.00 # line 2 length in m 
    >>> beta = 2*np.pi/3e8*frange # phase constant over freq
    >>> phase1 = -np.exp(1j*beta*line1_l)
    >>> phase2 = -np.exp(1j*beta*line2_l)
    >>> Smeas[:,0,0] = Smeas[:,0,0] * phase1 * phase1
    >>> Smeas[:,0,1] = Smeas[:,0,1] * phase1 * phase2
    >>> Smeas[:,1,0] = Smeas[:,1,0] * phase1 * phase2
    >>> Smeas[:,1,1] = Smeas[:,1,1] * phase2 * phase2
    >>> ### create an 2x2-Matrix N with gaussian noise #####
    >>> ns = len(frange) # number of samples
    >>> sigma = 0.01
    >>> np.random.seed(42)
    >>> N11 = np.random.normal(0,sigma,ns) + 1j * np.random.normal(0,sigma,ns) 
    >>> N12 = np.random.normal(0,sigma,ns) + 1j * np.random.normal(0,sigma,ns) 
    >>> N22 = np.random.normal(0,sigma,ns) + 1j * np.random.normal(0,sigma,ns) 
    >>> N = np.array([[N11,N12],[N12,N22]])
    >>> N = np.transpose(N,(2,0,1))
    >>> Smeas = Smeas + N
    >>> fig,ax1,ax2 = sparplot(frange,Smeas,True)
    >>> al = ax1.set_ylim(-30,2)
    >>> #ax2.plot(frange, np.angle(Smeas[:,1,1])*180/np.pi)
    >>> ### The filter is obviously detuned, now find the Coupling matrix using fit_to_M
    >>> ### use a butterworth 3-pole as initial guess
    >>> Mini = MofChebychev(n,0)
    >>> f0, BW, Qinitial = (2.0e9,0.1,1e10) ## De-Normalisations for CM-filter
    >>> ### Set a list of Coupling to be adjusted (optimized), 
    >>> ### here all direct couplings and self-couplings
    >>> Mindiceslist =[(1,0),(1,2),(2,3),(3,4),(1,1),(2,2),(3,3)]
    >>> ### since the response varies between 1.85 and 2.2 Hz only only 
    >>> ### we limit the optimisation to this range
    >>> fstart,fstop = (1.85e9,2.15e9)
    >>> ### to speed up the optimizer we can reduce the number 
    >>> ### of freq-samples by factor deci
    >>> deci = 2
    >>> Mfit, Qfit, phaseoffsetinput, phaseoffsetoutput, res = \
            fit_to_M(frange,Smeas,Mini,f0,BW,Qinitial,[],Mindiceslist, \
            False,deci,fstart,fstop,200) 
    >>> ### plot fitted response on top of measured response
    >>> Sfit = RespM2(Mfit, 2.0e9, 0.1, 1e9, frange)
    >>> Sfit[:,0,0] *=  np.exp(2j*(phaseoffsetinput*beta))
    >>> Sfit[:,0,1] *=  np.exp(1j*((phaseoffsetoutput+phaseoffsetinput)*beta))
    >>> fig,ax1,ax2 = sparplot(frange,Smeas,True)
    >>> fig.savefig("fitM1.png",dpi=300)
    >>> a = ax1.plot(frange, 20*np.log10(abs(Sfit[:,0,0])),'ro', \
            ms=2,label="S11 fitted")
    >>> a = ax1.plot(frange, 20*np.log10(abs(Sfit[:,0,1])),'ro',ms=2,label="S21 fitted")
    >>> a = ax1.legend()
    >>> a = ax1.set_ylim(-30,2)
    >>> a = ax2.plot(frange, np.angle(Sfit[:,0,0])*180/np.pi,'ro',ms=2,label="S11 fitted")
    >>> a = ax2.plot(frange, np.angle(Sfit[:,0,1])*180/np.pi,'ro',ms=2,label="S21 fitted")
    >>> a = ax2.legend()
    >>> fig.savefig("fitM2.png",dpi=300)
    >>> ### Compare Coupling Matrices 
    >>> print(np.around(Mfit,2))
    [[ 0.    0.6   0.    0.    0.  ]
     [ 0.6  -0.14  0.7   0.    0.  ]
     [ 0.    0.7   0.06  0.51  0.  ]
     [ 0.    0.    0.51 -0.21  0.7 ]
     [ 0.    0.    0.    0.7   0.  ]]
       
    .. figure:: fitM1.png
       :width: 100%
       
    .. figure:: fitM2.png
       :width: 100%
       
    .. raw:: latex
      
       \\clearpage
       \\newpage

    '''

    iter = 0
    n = np.shape(Minitial)[0]-2
    
    ##### Calculate deviation between measured S-par and fitted from M  S-par ########
    def differcalc(p,Mindiceslist,M,f0,BW,fdeci,s11deci,s21deci,s22deci):
        for kk in range(len(p)-3):
            i,j = Mindiceslist[kk]
            M[i,j] = p[kk]
            M[j,i] = p[kk]
            if is_symmetric: 
                M[n+1-i,n+1-j] =  p[kk]
                M[n+1-j,n+1-i] =  p[kk]
        if p[-3]<=0: # Q cannot be negative!
            Q = 10
        else:
            Q = p[-3]
        phaseoffset11 = p[-2]
        phaseoffset22 = p[-1]
        beta = 2*np.pi/3e8*fdeci
        S = RespM2(M, f0, BW, Q, fdeci, NRNlist)
        s11fit = S[:,0,0]
        s21fit = S[:,1,0]
        s22fit = S[:,1,1]
        s11fit *= np.exp(2j*phaseoffset11*beta)
        s21fit *= np.exp(1j*phaseoffset11*beta) * np.exp(1j*phaseoffset22*beta)
        s22fit *= np.exp(2j*phaseoffset22*beta)
        err =  np.concatenate ( ( \
               -np.log(abs(s11deci))*abs(s11deci-s11fit), \
               -np.log(abs(s21deci))*abs(s21deci-s21fit), \
               -np.log(abs(s22deci))*abs(s22deci-s22fit) ))
        #print(np.sqrt(np.sum(err**2)))
        return err
 
    #### Decimation to reduce Load #################################
    # findet index for sample at 10 MHz
    if fstart == 0:
        idx_start = 0
    else:
        idx_start  = (np.abs(flist-fstart)).argmin()
    if fstop == 0:
        idx_stop = len(flist)
    else:
        idx_stop  = (np.abs(flist-fstop)).argmin()
    
    #print(fstart,flist)
    s11meas,s12meas, s21meas, s22meas = (Smeas[:,0,0],Smeas[:,0,1],Smeas[:,1,0],Smeas[:,1,1]) 
    s11deci = s11meas[idx_start:idx_stop:deci]
    s21deci = s21meas[idx_start:idx_stop:deci]
    s22deci = s22meas[idx_start:idx_stop:deci]
    fdeci   =   flist[idx_start:idx_stop:deci]
    M = deepcopy(Minitial)
    
    plist = [M[ind] for ind in Mindiceslist]
    boundlower = [-10.]*len(plist)
    boundupper = [10.]*len(plist)
    
    ## Guess Phase Loading on boundaries #####
    #print("Phase Loading:")
    phase11 = (np.angle(s11meas))
    phase22 = (np.angle(s22meas))
    end = len(flist)
    ran = np.r_[0:10,end-11:end]
    pp11 = np.polyfit(flist[ran],phase11[ran],1)
    pp22 = np.polyfit(flist[ran],phase22[ran],1)
    phaseoffset11 = 3e8/2/np.pi*pp11[0] / 2.0 
    phaseoffset22 = 3e8/2/np.pi*pp22[0] / 2.0
    if phaseoffset11<0:phaseoffset11=0 
    if phaseoffset22<0:phaseoffset22=0 
    #print(phaseoffset11,phaseoffset22)
    plist.append(Qinitial)
    plist.append(phaseoffset11)
    plist.append(phaseoffset22)
    boundlower = boundlower + [0,-500.,-500.]
    boundupper = boundupper + [1e99,500.,500.]
    res = least_squares(differcalc, plist,xtol = 1e-12,loss='soft_l1', bounds=(boundlower,boundupper),args=(Mindiceslist,M,f0,BW,fdeci,s11deci,s21deci,s22deci)) 
    # print("** Done** - Cost of Fitting =",res.cost)
    if res.cost>maxerr/deci:
        print(res.message)
        raise RuntimeError("M-Fit did not succeed, no Fitting found, Cost="+str(res.cost)) 
    return M,res.x[-3],res.x[-2],res.x[-1],res 


########################################################################
def showmatrixplot(M,color=True):
    '''
    Creates a colorful plot of the Coupling Matrix
    
    Parameters
    ----------
    M : 2-dim array of floats
        Coupling Matrix n+2 x n+2

   
    Examples:
    ---------
    >>> M = MofChebychev(5,0.1)
    >>> M = np.array([[0.   , 0.6  , 0.   , 0.   , 0.   ], \
                   [0.6  , -0.14 , 0.705, 0.   , 0.   ], \
                   [0.   , 0.705, 0.06 , 0.505, 0.   ], \
                   [0.   , 0.   , 0.505, -0.211 , 0.703], \
                   [0.   , 0.   , 0.   , 0.703, 0.   ]])
    >>> M[2,2] = -0.2
    >>> M[1,4] = -0.7
    >>> fig, ax = showmatrixplot(M)
    >>> fig.savefig('matrixplot.png',dpi=300)
    
    .. figure:: matrixplot.png
       :width: 50%
 
    '''
    
    if color:
        vv = 0.5
    else:
        vv = 500
    fig,ax = plt.subplots(figsize=(6,6))
    res = ax.imshow((M), vmin=-vv, vmax=vv, cmap=plt.cm.seismic, 
                interpolation='nearest')
    width, height = M.shape
    for x in range(width):
      for y in range(height):
        ax.annotate(str(np.around(M[x,y],3)), xy=(y, x), 
                    horizontalalignment='center',
                    verticalalignment='center',fontsize=14)
    return fig, ax

        

#######################################################################
def sparplot(frange,S,phase = False):
    '''
    Create an S-parameter freq.-response plot of an f/S pair list
    
    Examples:
    ---------
    >>> ripple = 1.0
    >>> M = MofChebychev(5,ripple)
    >>> frange = np.arange(1.8,2.3,0.001)
    >>> S = RespM2(M, 2.0, 0.1, 1e9, frange)
    >>> fig, ax = sparplot(frange,S)
    >>> ll = ax.set_ylim(-20,1)
    >>> lh = ax.axhline(-ripple,color='r',linestyle=":");
    >>> lt = ax.text(2.15,-ripple-1.1,"Ripple "+str(ripple)+"dB",fontsize=14,color='r');
    >>> fig.savefig('sparplot.png',dpi=300)
    
    .. figure:: sparplot.png
       :width: 80%
        
    .. raw:: latex
      
       \\clearpage
       \\newpage

    '''
    ### in dB #############
    s11 = S[:,0,0]
    s21 = S[:,1,0]
    s11dB = [20*np.log10(abs(x)) for x in s11]
    s21dB = [20*np.log10(abs(x)) for x in s21]
    s11phase = (np.angle(s11))*180/np.pi 
    s21phase = (np.angle(s21))*180/np.pi
    
    if phase:
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4))
        ax1.plot(frange,s11dB, label="$S_{11}$")
        ax1.plot(frange,s21dB, label="$S_{21}$")
        ax1.set_ylim(-40.0, 2.0)
        ax1.set_xlabel("f in Hz")
        ax1.set_ylabel("S in dB")
        ax1.grid()
        ax2.set_title("Magnitude")
        ax1.legend()
        ax2.plot(frange,s11phase, label="$S_{11}$")
        ax2.plot(frange,s21phase, label="$S_{21}$")
        #ax2.set_ylim(-180.0, 180.0)
        ax2.set_xlabel("f in Hz")
        ax2.set_ylabel("S in degrees")
        ax2.grid()
        ax2.set_title("Phase")
        ax2.legend()
        return fig,ax1,ax2
    else:
        fig,ax = plt.subplots(figsize=(6,4))
        ax.plot(frange,s11dB, label="$S_{11}$")
        ax.plot(frange,s21dB, label="$S_{21}$")
        ax.set_ylim(-40.0, 2.0)
        ax.set_xlabel("f in Hz")
        ax.set_ylabel("S in dB")
        ax.grid()
        ax.legend()
        return fig,ax


###############################################################################    
def drawcouplingdiagram(M,respos,ntype): 
    '''
    Draw the coupling diagram in bullet nodes and line link style
    
    Parameters
    ----------
   
       M: 2dim-square-array
            coupling matrix 
       respos: (int,int)
          positions of the resonators as 2dim tuples
       ntype: string
          resonator type, Port 'P', resonaotr 'R' or NR-Node 'NRN'
       
    Returns
    -------
    
    drawing object
         d
        

    Examples
    --------

    >>> M = MofChebychev(4,0.5)
    >>> M[1,4] = 0.111
    >>> M[4,4] = -0.1
    >>> respos = [(0,0),(1,0),(2,0),(2,1),(1,1),(0,1)]
    >>> ntype =  ['P',  'R',  'R',  'R',  'R',  'P']
    >>> d = drawcouplingdiagram(M,respos,ntype)
    >>> d.draw(showplot=False) # Use showplot = True for interactive Plots
    >>> d.save('couplingdiagram.png',dpi=300)
    >>> plt.close(d.fig)

    .. image:: couplingdiagram.png
    
    '''
    
    # Resoposition and res type ntype: R,P, or NRN
    respos = np.array(respos)
    n = np.size(M,0)-2
    if n+2>len(respos):
        print("M-Order+2 higher than number of placed Resonators",n+2,len(respos))
        return
    if len(ntype) !=len(respos):
        print("size ntype does not match size respos",len(ntype),len(respos))
        return
    
    d = schem.Drawing()
    for pp in range(len(respos)):
        xy = respos[pp]*3
        if ntype[pp]=='P':
            d.add(PORT,xy=xy)
        elif ntype[pp]=='R':
            d.add(RESO,xy=xy)
        elif ntype[pp]=='NRN':
            d.add(NRN,xy=xy)

    for i in range(n+2):      ### Make symmetric
        for j in range(i,n+2):
            M[j,i] = M[i,j]
            if M[i,j]: 
                m = d.add(e.LINE,xy=respos[i]*3,to=respos[j]*3)
                m.add_label(str(round(M[i,j],3)),size=14)
    return d



if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
 

