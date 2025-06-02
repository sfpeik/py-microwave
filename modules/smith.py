# Plot a smith chart in a given (ax) plot
###########################################
"""This module allows the generation of Smith charts"""

__license__ = "Soeren Peik"
__revision__ = " 2021-April "
__docformat__ = 'reStructuredText'


from itertools import count

import matplotlib.pyplot as plt
from matplotlib import patches
from numpy import (real, imag, logspace, linspace, log10, exp, pi, angle, arange, tan, append,
                   cos, arccos, arcsin, arctan2, sin, around, array, nan, sqrt, conj, argmin, concatenate)

import schemdraw as schem
import schemdraw.elements as e

import os
import svgutils.compose as sc

#from yt_dlp.utils import find_xpath_attr

### convert mag, phase into complex number
### phase is given in degrees

def magphase(amp, phase):
    return amp * exp(1j * phase / 180.0 * pi)


def magphase_str(c):
    """
    Returns a nicely formatted String to print complex numbers in ampl. and phase
    uses Unicode characters!
    """
    return u'{0:6.3f}\u2220{1:5.1f}\u00B0'.format(abs(c), angle(c) * 180 / pi)


def magphase_tex(c):
    """
    Returns a TeX formatted String to print complex numbers in ampl. and phase
    """
    return r' {0:6.3f} \angle {1:5.1f}^\circ '.format(abs(c), angle(c) * 180 / pi)


def magphasetex(c):  # alternative without underscore
    return magphase_tex(c)


def reflcoeff(Z, Z0=50):
    r"""
    Return the reflection coefficient for a given load impedance Z and line impedance Z0

    Returns:

    .. math::
        \Gamma = \\frac{Z-Z_0}{Z+Z_0}

    :example:

        .. code-block:: python

          Z = 100.
          Z0 = 50.
          gam = smith.ReflCoeff(Z,Z0)
          print(gam)

    """
    return (Z - Z0) / (Z + Z0)

hide = False
def hidesolution(status= True):
    hide = status

### convert mag, phase into complex number
### phase is given in degrees

########################################################################
########################################################################
class Element:
    _ids = count(0)

    def __init__(self, elementtype, connection, value, Zin=nan, attributes=None):
        if attributes is None:
            attributes = {}
        self.type = elementtype
        self.connection = connection
        self.value = value
        self.atFreq = None
        self.compValue = None
        self.id = next(self._ids)
        self.attr = attributes
        self.Zin = Zin

    def __str__(self):
        line = "Element:{0:2} {1:12} {2:10} Z={3:10} {4:30}".format(self.id, self.type, self.connection, self.value,
                                                                    str(self.attr))

        return line


########################################################################
########################################################################
def gam(Z, Z0):
    return (Z - Z0) / (Z + Z0)


# Superimpose a paper chart onto the plot ##############################
def superimposeChart(fig, svgfile='smith_paper.svg'):
    fig.set_size_inches(18, 18)
    plt.tight_layout(pad=0.0)
    fig.savefig('smithdata.svg', transparent=True)
    sc.Figure("350mm", "350mm",
              sc.Panel(sc.SVG(os.path.dirname(__file__)+"/smith_mwl_chartonly.svg").scale(1.27).move(174, 195)),
              sc.Panel(sc.SVG("smithdata.svg").scale(1.0).move(0, 0))
              ).save(svgfile)
    print("Saved as ",svgfile)
    return svgfile

########################################################################################################################
########################################################################################################################


class Smith:
    def __init__(self, ax, typ='smith', Z0=1, fontscale=None, trans= None, fineness=1, **kwargs):
        
        self.scaleadjust = ax.get_window_extent().height / 1000.
        if fontscale is None:
            self.fontscale = 1* sqrt(self.scaleadjust)
        else:
            if fontscale >10:
                raise ValueError("Deprecated fontscale in use, must be smaller than 10.0, here: "+str(fontscale))
            self.fontscale = fontscale* sqrt(self.scaleadjust)
        self.linefactor = 1.0
        if 'lw' in kwargs:
            self.linefactor = kwargs['lw']
            del kwargs['lw']
        if 'linewidth' in kwargs:
            self.linefactor = kwargs['linewidth']
            del kwargs['linewidth']
        ### Set Transparency legacy style (deprecated)
        if not 'alpha' in kwargs:
            kwargs['alpha'] = 0.5
        elif trans is not None:
            kwargs['alpha'] = trans
        self.linewidth = 2

        # print("Scale Factor:", self.scaleadjust)
        Element._ids = count(1)
        self.clist = [[], [], [], []]
        self.cendlist = [[], [], [], []]
        self.cshow = [[], [], [], []]
        self.cshowx = [[], [], [], []]
        self.clist[0] = [0, 0.5, 1, 2]
        self.cendlist[0] = [100, 100, 100, 100]
        self.cshow[0] = [1, 1, 1, 1]
        self.cshowx[0] = [0, 1, 1, 1]

        self.clist[1] = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 4, 10, 20]
        self.cendlist[1] = [100, 2, 1, 2, 1, 100, 2, 10, 10, 10, 100]
        self.cshow[1] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        self.cshowx[1] = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        self.clist[2] = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0,
                         5.0, 6, 7, 8, 9, 10, 15, 20, 50]
        self.cendlist[2] = [100, 2, 2, 2, 2, 5, 2, 2, 2, 2, 100, 2, 2, 2, 2, 50, 5, 5, 5, 5, 5, 5, 5, 5, 50, 50, 50, 50]
        self.cshow[2] = [0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1]
        self.cshowx[2] = [0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1]

        self.clist[3] = self.clist[2] + list(arange(0.05, 0.5, 0.1))
        self.cendlist[3] = self.cendlist[2] + [0.5] * 5
        self.cshow[3] = self.cshow[2] + [0] * 5

        self.clist = self.clist[fineness]
        self.cendlist = self.cendlist[fineness]
        self.cshow = self.cshow[fineness]

        self.cfat = [0, 0.2, 0.5, 1, 2, 5, 10]
        self.f = 0
        self.ax = ax
        self.Z0 = Z0
        self.kwargs = kwargs
        self.elements = []
        self.impedances = []
        self.elementvalues = []
        self.elementvaluesunit = []
        self.circuittypes = []
        self.circuit = []
        self.Zinlist = []
        self.isnorm = False
        self.paper = False
        self.showArrows = False
        e.style(e.STYLE_IEC)

        if Z0 == 1:
            self.isnorm = True
        if typ == 'smith' or type == 'impedance':
            if not ('c' in self.kwargs or 'color' in self.kwargs):
                self.kwargs['color'] = 'r'
            self.addimpedancegrid()
        elif typ == 'inverted' or typ == 'admittance':
            if not ('c' in self.kwargs or 'color' in self.kwargs):
                self.kwargs['color'] = 'b'
            self.addadmittancegrid()
        elif typ == 'both':
            if not ('c' in self.kwargs or 'color' in self.kwargs):
                self.kwargs['color'] = 'r'
            self.addimpedancegrid()
            if self.kwargs['color']=='r':
                self.kwargs['color'] = 'b'
            self.addadmittancegrid()
        elif typ == 'paper':
            self.addemptygrid()
            self.paper = True
        else:
            raise ValueError('Smith Style not defined')

    def addemptygrid(self,visible = False):
        self.ax.set_aspect('equal', anchor='C')
        self.fontscale = self.fontscale *2
        self.ax.plot(sin(linspace(0, 2. * pi)), cos(linspace(0, 2. * pi)), 'k--', lw=1.) # plot outline
        self.ax.plot(1.5*sin(linspace(0, 2. * pi)), 1.5*cos(linspace(0, 2. * pi)), 'k--', lw=0.)  # plot outline invisible
        if not visible:
            self.ax.axis('off')

    def addimpedancegrid(self, inverted = False, **kwargs):
        '''
        :param inverted: plots an inverted chart, aka Admittance Chart
        :param kwargs: plot arguments, e.g. color, lw, ...
        :return:
        '''
        if not kwargs:
            kwargs = self.kwargs
        con = -1
        color = "r"
        if 'color' in self.kwargs:  color = self.kwargs['color']
        else: self.kwargs['color'] = color
        if 'alpha' in self.kwargs: alpha = self.kwargs['alpha']

        ## Cirlces ###########
        for zreal in self.clist:
            con = con + 1
            zimag = logspace(-4, log10(500), 200)
            z = zreal + 1j * zimag
            gamma = (z - 1) / (z + 1)
            if inverted:
                gamma = - gamma
            lab = ' {:1.1f}'.format(zreal).replace('.0', '')
            ### plots ################
            if zreal in self.cfat:
                llw = 1.0 * self.linefactor
            else:
                llw = 0.4 * self.linefactor
            self.ax.plot(real(gamma), imag(gamma), linewidth=llw, **kwargs)
            self.ax.plot(real(gamma), -imag(gamma), linewidth=llw, **kwargs)

            ### Labels ##################
            if self.cshow[con]:
                gamlabel = (zreal - 1.0) / (zreal + 1.0)
                orient = 1
                if inverted:
                    gamlabel = -gamlabel
                    orient = -1
                self.ax.text(real(gamlabel) , 0.0, lab,rotation_mode='anchor', ha='left', va='bottom',
                             rotation=90*orient, fontsize=10 * self.fontscale, color=color, alpha = alpha)
            if zreal in self.cfat:
                gamlabel = (zreal + 1j - 1) / (zreal + 1j + 1)
                anglabel = 180 + angle(gamlabel - 1 - 1j) * 180 / pi
                if inverted: gamlabel = -gamlabel
                if zreal != 0:
                    self.ax.text(real(gamlabel), imag(gamlabel), lab,rotation_mode='anchor',  ha='left', va='bottom',
                                 rotation=anglabel, fontsize=10 * self.fontscale,  color = color, alpha = alpha)
                    self.ax.text(real(gamlabel), imag(-gamlabel), lab,rotation_mode='anchor', ha='left', va='bottom',
                                 rotation=180-anglabel, fontsize=10 * self.fontscale, color = color, alpha = alpha)
        con = -1
        ## Bends #############  
        for zimag in self.clist:
            con = con + 1
            zreal = logspace(-4, log10(self.cendlist[con]), 100)
            z = zreal + 1j * zimag
            gamma = (z - 1) / (z + 1)
            if inverted: gamma = -gamma
            if zimag in self.cfat:
                llw = 0.8 * self.linefactor
            else:
                llw = 0.4 * self.linefactor
            self.ax.plot(real(gamma), imag(gamma),  linewidth=llw, **kwargs)
            self.ax.plot(real(gamma), -imag(gamma),  linewidth=llw, **kwargs)
            ### add labels ###############
            if zimag == 0: continue
            lab = ' %3.1fj ' % zimag
            gg = (1j * zimag - 1) / (1j * zimag + 1)
            if inverted: gg = -gg
            x = real(gg)
            y = imag(gg)
            ang = angle(gg) * 180.0 / pi
            if self.cshow[con]:
                self.ax.text(x, y, lab,  ha='right', va='bottom', rotation_mode='anchor', rotation=ang,
                             fontsize=10 * self.fontscale, alpha = alpha, color = color)
            lab = ' -%3.1fj' % zimag
            if self.cshow[con]:
                self.ax.text(x, -y, lab,  ha='left', va='bottom',rotation_mode='anchor', rotation=180 - ang,
                             fontsize=10 * self.fontscale, alpha = alpha, color = color)
        self.ax.set_xlim(-1.2, 1.2)
        self.ax.set_ylim(-1.2, 1.2)
        self.ax.axis('off')
        self.ax.axis('equal')

    ##### plot an inverted chart in blue

    def addadmittancegrid(self):
        self.addimpedancegrid(inverted = True)

    def addpolargrid(self, **kwargs):
        '''
        Add a polar grid around the edge of the chart
        :param kwargs:
        :return:
        '''

        labelradius = 1.2
        if not "color" in kwargs:   kwargs['color'] = 'k'
        if not "lw"    in kwargs:  kwargs['lw'] = 0.5
        if not "alpha" in kwargs:  kwargs['alpha'] = 0.5
        for magn in arange(0.1, 1, 0.1):
            patch = plt.Circle((0, 0), magn, fill=False, **kwargs)
            self.ax.add_patch(patch)
            lab = "{:1.1f}".format(magn)
            self.ax.text(0, magn + 0.05, lab, ha='center', va='top',
                         rotation=0, fontsize=10 * self.fontscale, color = kwargs['color'], alpha = kwargs['alpha'])
        for ang in arange(0, 360, 10):
            xl = cos(ang * pi / 180)
            yl = sin(ang * pi / 180)
            self.ax.plot((0, labelradius * xl), (0, labelradius * yl),  **kwargs)
            lab = "{:3d}°".format(ang)
            if ang > 180: lab = "{:3d}°".format(-360 + ang)
            xl = labelradius * cos(ang * pi / 180 + 0.025)
            yl = labelradius * sin(ang * pi / 180 + 0.025)
            self.ax.text(xl, yl, lab, ha='center', va='center', rotation=ang,
                         fontsize=10 * self.fontscale, color = kwargs['color'], alpha = kwargs['alpha'])

    ##### Add a ring of angle Values #################################################################

    def addanglering(self, **kwargs):
        '''
        Adds a ring with angles in lambdas to chart
        :param color:
        :param args:
        :param kwargs:
        :return:
        '''
        if not "color" in kwargs: kwargs['color'] = 'k'
        if not "lw"    in kwargs: kwargs['lw'] = 1
        if not "alpha" in kwargs: kwargs['alpha'] = 0.5

        for lam in arange(0.0, 0.5, 0.02):
            ang = pi - 4 * pi * lam
            x = 1.1 * cos(ang)
            y = 1.1 * sin(ang)
            self.ax.plot((0.95 * x, 1.05 * x), (0.95 * y, 1.05 * y), **kwargs)
            lab = str(lam) + r'$\lambda$'
            self.ax.text(x*0.985 , y*0.985 , lab, ha='left', va='bottom',rotation_mode='anchor', rotation=ang * 180.0 / pi,
                         fontsize=10 * self.fontscale, color = kwargs['color'], alpha = kwargs['alpha']   )
        for lam in arange(0.0, 0.5, 0.002):  # minor tics
            ang = pi - 4 * pi * lam
            x = 1.1 * cos(ang)
            y = 1.1 * sin(ang)
            length1 = 1.00
            length2 = 0.97
            if int(lam * 1000) % 10 == 0:
                self.ax.plot((0.95 * x, length1 * x), (0.95 * y, length1 * y), **kwargs)
            else:
                self.ax.plot((0.95 * x, length2 * x), (0.95 * y, length2 * y), **kwargs)
        #self.ax.axis('off')
        self.ax.axis('equal')
        plt.tight_layout(pad=.0)

    ##### Add a Ruler ##################################################################################

    def addruler(self, *args, **kwargs):
        '''
        Adds a ruler underneath the chart
        :param args:
        :param kwargs:
        :return:
        '''
        mfs = 7
        #### Gamma Values #####
        yl = -1.6
        self.ax.text(0.04, yl - 0.05, "Refl Coeff $\\Gamma$", color='k', ha='left', va='center',
                     fontsize=10 * self.fontscale, *args, **kwargs)
        self.ax.plot((-1, 0), (yl, yl), 'k', linewidth=1, *args, **kwargs)
        for xl in arange(0, 1.1, 0.1):
            self.ax.plot((-xl, -xl), (yl, yl - 0.02), 'k', linewidth=1, *args, **kwargs)
            lab = "{:1.1f}".format(xl)
            self.ax.text(-xl, yl - 0.05, lab, color='k', ha='center', va='center', fontsize=mfs * self.fontscale, *args,
                         **kwargs)
        for xl in arange(0, 1.0, 0.02):
            self.ax.plot((-xl, -xl), (yl, yl - 0.01), 'k', linewidth=0.5, *args, **kwargs)

        #### dB Values ######
        yl = -1.5
        self.ax.text(0.04, yl - 0.05, "Return Loss RL", color='k', ha='left', va='center',
                     fontsize=10 * self.fontscale)
        self.ax.plot((-1, 0), (yl, yl), 'k', linewidth=1)
        for dB in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 1e20]:
            xl = 10 ** (-dB / 20)
            self.ax.plot((-xl, -xl), (yl, yl - 0.02), 'k', linewidth=1)
            lab = "{:1.0f}".format(dB)
            if dB > 1e10: lab = "$\\infty$"
            self.ax.text(-xl, yl - 0.05, lab, color='k', ha='center', va='center', fontsize=mfs * self.fontscale)

        #### VSWR Values ######
        yl = -1.4
        self.ax.text(0.04, yl - 0.05, "VSWR", color='k', ha='left', va='center', fontsize=10 * self.fontscale)
        self.ax.plot((-1, 0), (yl, yl), 'k', linewidth=1)
        for VSWR in [1, 1.2, 1.5, 1.8, 2, 2.5, 3, 4, 5, 10, 20, 1e20]:
            xl = (VSWR - 1) / (VSWR + 1)
            self.ax.plot((-xl, -xl), (yl, yl - 0.02), 'k', linewidth=1)
            lab = "{:1.1f}".format(VSWR)
            if VSWR > 1e10: lab = r"$\infty$"
            self.ax.text(-xl, yl - 0.05, lab, color='k', ha='center', va='center', fontsize=mfs * self.fontscale)

        plt.tight_layout(pad=.0)

    ##### Add a ruler marker ################################################################
    def addrulermarker(self, mag, color='b', *args, **kwargs):
        '''
        Add a marker across the lower ruler.
        :param mag: the magnitude of the Gamma value to be marked
        :param color:
        :param args:
        :param kwargs:
        :return:
        '''
        xl = abs(mag)
        if self.paper:
            yl = -1.45
            yl2 = -1.2
        else:
            yl = -1.75
            yl2 = -1.45
            plt.tight_layout(pad=2.0)
        self.ax.plot((-xl, -xl), (yl, yl2), color=color, linewidth=1.6, *args, **kwargs)


    ##### Add a ruler distance ################################################################
    def addrulerdistance(self, mag, *args, **kwargs):
        yl = -1.6
        xl = abs(mag)
        self.ax.plot((0, -xl), (yl, yl), 'b', linewidth=3, *args, **kwargs)
        self.ax.plot((-xl, -xl), (yl - 0.02, yl + 0.02), 'b', linewidth=3, *args, **kwargs)
        self.ax.plot((0, 0), (yl - 0.02, yl + 0.02), 'b', linewidth=3, *args, **kwargs)
        plt.tight_layout(pad=2.0)

    ##### Add a circle in a smith chart ################################################################

    def addcircle(self, C, R, color='y', alpha=0.3, **args):
        '''
        add a circle to the chart
        :param C:  center of circle
        :param R: radius of circle
        :param color: color
        :param alpha: transparency
        :param args:
        :return:
        '''
        patch = plt.Circle((real(C), imag(C)), R, alpha=1.0, fill=False, ec=color, linewidth=2, **args)
        self.ax.add_patch(patch)
        patch = plt.Circle((real(C), imag(C)), R, alpha=alpha, fc=color, linewidth=0)
        self.ax.add_patch(patch)

    ##### Add an arrow in a smith chart ################################################################

    def addarrow(self, c, color='k', *args, **kwargs):
        ang = angle(c)
        x = real(c) - 0.1 * cos(ang)
        y = imag(c) - 0.1 * sin(ang)
        self.ax.arrow(0, 0, x, y, head_width=0.03, head_length=0.1, alpha=1, fc=color, ec=color, *args, **kwargs)

    ##### Add a point in a smith chart ###############################################################

    def addgammapoint(self, Gamma, label='', ori='NE', **kwargs):
        '''
        Adds a marker at the given reflection coefficient
        the label can be set
        :param Z:
        :param label: Label for the Marker
        :param ori: Location of the Label in compass direction, NE,NW,SE, or SW
        :param kwargs: Setting other attributes of point
        :return:
        '''
        g = Gamma
        lab = label
        if not "color" in kwargs:  kwargs['color'] = 'k'
        if not "lw" in kwargs:  kwargs['lw'] = 1
        if not "alpha" in kwargs:  kwargs['alpha'] = 0.9
        patch = plt.Circle((real(g), imag(g)), 0.012, **kwargs)
        self.ax.add_patch(patch)
        orix, oriy = (0, 0)
        if ori == 'NE': orix = 1; oriy = 1
        if ori == 'NW': orix = -2.3; oriy = 1
        if ori == 'SE': orix = 1; oriy = -1.6
        if ori == 'SW': orix = -2.3; oriy = -1.6
        self.ax.annotate(lab, xy=(real(g), imag(g)), xycoords='data', fontsize=12 * self.fontscale,
                         xytext=(30 * orix * self.fontscale, 30 * oriy * self.fontscale), textcoords='offset points',
                         bbox=dict(boxstyle="round", fc="0.8", alpha=1.0),
                         arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=90,rad=20"), )


    def addpoint(self, Z, label='Z', ori='NE', **kwargs):
        '''
        Adds a marker at the given impedance
        the label includes the Z and Y Value at the point normalized
        :param Z:
        :param label: Label for the Impedance
        :param ori: Location of the Label in compass direction, NE,NW,SE, or SW
        :param kwargs: Setting other attributes of point
        :return:
        '''
        g = gam(Z, self.Z0)
        if not "color" in kwargs:  kwargs['color'] = 'k'
        if not "lw" in kwargs:  kwargs['lw'] = 1
        if not "alpha" in kwargs:  kwargs['alpha'] = 0.9
        patch = plt.Circle((real(g), imag(g)), 0.012, **kwargs)
        self.ax.add_patch(patch)
        Z = around(real(Z), 3) + 1j * around(imag(Z), 3)  # fix small non-zero values
        if abs(Z) < 1e-10:  # short
            lab = '%s = 0 \n $Y$ = $\\infty$' % label
        elif abs(Z) > 1e20:  # open
            lab = '%s = $\\infty$ \n $Y$ = 0 ' % label
        else:
            Y = 1 / Z
            lab = '%s = %4.2f%+4.2fj \n $Y$ =  %4.2g%+4.2fj' % (
                label, real(Z / self.Z0), imag(Z / self.Z0), real(self.Z0 * Y), imag(self.Z0 * Y))
        orix, oriy = (0, 0)
        if ori == 'NE': orix = 1; oriy = 1
        if ori == 'NW': orix = -4.3; oriy = 1
        if ori == 'SE': orix = 1; oriy = -1.6
        if ori == 'SW': orix = -4.3; oriy = -1.6
        self.ax.annotate(lab, xy=(real(g), imag(g)), xycoords='data', fontsize=12 * self.fontscale,
                         xytext=(30 * orix * self.fontscale, 30 * oriy * self.fontscale), textcoords='offset points',
                         bbox=dict(boxstyle="round", fc="0.8", alpha=1.0),
                         arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=90,rad=20"), )


    ##### Add a Starting Point    ###############################################################
    def addstart(self, Z1, **kwargs):
        if not "color" in kwargs: kwargs['color'] = 'k'
        if not "lw" in kwargs: kwargs['lw'] = 3
        if not "alpha" in kwargs: kwargs['alpha'] = 0.9
        Znew = Z1
        g = gam(Z1, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, color = kwargs['color'])
        self.ax.add_patch(patch)
        self.appendelement('Start', Z1, Znew)
        return Znew

    def addinput(self, label=""):
        self.appendelement('inputport',1,1)
        return

    ##### Add impedance in series ###############################################################

    def addseries(self, Z1, Z2, **kwargs):
        if not "color" in kwargs: kwargs['color'] = 'r'
        if not "lw" in kwargs: kwargs['lw'] = self.linewidth
        Znew = Z1 + Z2
        t = arange(0.01, 1, 0.01)
        Z = Z1 + Z2 * t
        self.ax.plot(real(gam(Z, self.Z0)), imag(gam(Z, self.Z0)), **kwargs)
        g = gam(Znew, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, alpha=0.9, color = kwargs['color'],lw = 0)
        self.ax.add_patch(patch)
        self.appendelement('Series', Z2, Znew)
        return Znew

    ##### Add impedance in parallel ###############################################################

    def addpara(self, Z1, Z2, **kwargs):
        if not "color" in kwargs: kwargs['color'] = 'b'
        if not "lw" in kwargs: kwargs['lw'] = self.linewidth
        Znew = Z1 * Z2 / (Z1 + Z2)
        t = arange(0.01, 1, 0.01)
        Z = (Z1 * Z2 / t) / (Z1 + Z2 / t)
        self.ax.plot(real(gam(Z, self.Z0)), imag(gam(Z, self.Z0)), **kwargs)
        g = gam(Znew, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, alpha=0.9, color = kwargs['color'], lw = 0 )
        self.ax.add_patch(patch)
        self.appendelement('Parallel', Z2, Znew)
        return Znew

    ##### Append type for lists  ###############################################################

    def appendelement(self, connection, Z, Znew, attributes=None):
        if attributes is None:
            attributes = {}
        epsilon = 1e-9
        if 'insert' in connection:
            typus = 'Line'
        elif 'stubopen' in connection:
            typus = 'Open Stub'
            connection = 'parallel'
        elif 'stubshort' in connection:
            typus = 'Shorted Stub'
            connection = 'parallel'
        elif 'port' in connection:
            typus = 'Port'
        elif abs(Z) < epsilon:
            typus = 'Short'
        elif abs(Z) > 1 / epsilon:
            typus = 'Open'
        elif real(Z) == 0 and imag(Z) < 0:
            typus = 'Capacitor'
        elif real(Z) == 0 and imag(Z) > 0:
            typus = 'Inductor'
        elif real(Z) > 0 and imag(Z) == 0:
            typus = 'Resistor'
        else:
            typus = 'Impedance'
        newelement = Element(typus, connection, Z, Znew, attributes)
        self.circuit.append(newelement)
        #print(newelement)

    ##### Insert Line ###############################################################

    def addline(self, Z1, length, Zline=50, **kwargs):
        if not "color" in kwargs: kwargs['color'] = 'g'
        if not "lw" in kwargs: kwargs['lw'] = self.linewidth

        Znew = Zline * (Z1 + 1j * Zline * tan(2 * pi * length)) / (Zline + 1j * Z1 * tan(2 * pi * length))
        t = arange(0.005, 1, 0.005)
        Z = Zline * (Z1 + 1j * Zline * tan(2 * pi * t * length)) / (Zline + 1j * Z1 * tan(2 * pi * t * length))
        self.ax.plot(real(gam(Z, self.Z0)), imag(gam(Z, self.Z0)), **kwargs)
        g = gam(Znew, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, alpha=0.9,color = kwargs['color'], lw = 0)
        self.ax.add_patch(patch)
        self.appendelement('inserted', Zline, Znew, {'length': length})
        return Znew

    ##### Draw Angle marker ###############################################################

    def addangle(self, Zx):
        # angle=pi-lam*4*pi
        ang = angle(gam(Zx, self.Z0))
        lam = round((pi - ang) / 4 / pi, 3)
        x = 1.2 * cos(ang)
        y = 1.2 * sin(ang)
        xl = 1.25 * cos(ang + 0.025 * self.fontscale)
        yl = 1.25 * sin(ang + 0.025 * self.fontscale)
        self.ax.plot((0, x), (0, y), 'k', linewidth=1)
        lab = str(lam) + r'$\lambda$'
        bbox_props = dict(boxstyle="round4", fc=(0.95, 0.95, 0.8), ec="0.1", alpha=1)
        self.ax.text(xl, yl, lab, color='k', ha='center', va='center', rotation=ang * 180.0 / pi,
                     bbox=bbox_props, fontsize=12 * self.fontscale)
        self.ax.set_xlim(-1.0, 1.0)
        self.ax.set_ylim(-1.0, 1.0)
        self.ax.axis('equal')
        plt.tight_layout(pad=2.0)

    def drawangle(self, Zx):
        print("drawangle is deprecated, use addangle ")
        self.addangle(Zx)

    ##### Add Open Stub ###############################################################

    def addstubopen(self, Z1, length, Zline):
        Zstub = -1j * Zline / tan(2 * pi * length)
        Znew = Z1 * Zstub / (Z1 + Zstub)
        t = arange(0.01, 1, 0.002)
        Zstub = -1j * Zline / tan(2 * pi * t * length)
        Z = Z1 * Zstub / (Z1 + Zstub)
        self.ax.plot(real(gam(Z, self.Z0)), imag(gam(Z, self.Z0)), 'y', linewidth=2)
        g = gam(Znew, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, alpha=0.9, linewidth=1)
        self.ax.add_patch(patch)
        self.appendelement('stubopen', Zline, Znew, {'length': length})
        return Znew

    ##### Add short Stub ###############################################################

    def addstubshort(self, Z1, length, Zline):
        Zstub = 1j * Zline * tan(2 * pi * length)
        Znew = Z1 * Zstub / (Z1 + Zstub)
        t = arange(0.01, 1, 0.01)
        Zstub = 1j * Zline * tan(2 * pi * length) / t
        Z = Z1 * Zstub / (Z1 + Zstub)
        self.ax.plot(real(gam(Z, self.Z0)), imag(gam(Z, self.Z0)), 'y', linewidth=2)
        g = gam(Znew, self.Z0)
        patch = plt.Circle((real(g), imag(g)), 0.01, alpha=0.9, linewidth=1)
        self.ax.add_patch(patch)
        self.appendelement('stubshort', Zline, Znew, {'length': length})
        return Znew

    ######## Return the schematic as plot object #########################################
    def plotschematic(self, fontsize = 12):
        if self.isnorm:
            unit = ''
        else:
            unit = r'\Omega'
            # plt.ioff()
        d = schem.Drawing(fontsize=fontsize)

        for ele in self.circuit:
            connection = ele.connection
            typus = ele.type
            impedance = ele.value
            if connection == 'Start':
                if typus == "Capacitor":
                    startC(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
                elif typus == "Inductor":
                    startL(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
                else:
                    startR(d, '$' + str(impedance) + unit + '$')
            if connection == 'inputport':
                inputport(d, "Zin\n\n\n\n")
                return d
            if typus == "Capacitor" and connection == "Parallel":
                shuntC(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
            elif typus == "Inductor" and connection == "Parallel":
                shuntL(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
            elif typus == "Resistor" and connection == "Parallel":
                shuntR(d, '$' + str(impedance) + unit + '$')
            elif typus == "Impedance" and connection == "Parallel":
                shuntR(d, '$' + str(impedance) + unit + '$')    
            elif typus == "Capacitor" and connection == "Series":
                seriesC(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
            elif typus == "Inductor" and connection == "Series":
                seriesL(d, '$' + str(imag(impedance)) + 'j' + unit + '$')
            elif typus == "Resistor" and connection == "Series":
                seriesR(d, '$' + str(impedance) + unit + '$')
            elif typus == "Impedance" and connection == "Series":
                seriesR(d, '$' + str(impedance) + unit + '$')    

            if typus == "Line":
                lamlen = ele.attr['length']
                impedance = ele.value
                lamstr = str(round(ele.attr['length'],3)) + r'\lambda'
                if abs(lamlen - 1. / 8) < 0.001:
                    lamstr = r'\lambda/8'
                if abs(lamlen - 1. / 4) < 0.001:
                    lamstr = r'\lambda/4'
                if abs(lamlen - 1. / 3) < 0.001:
                    lamstr = r'\lambda/3'
                if abs(lamlen - 1. / 6) < 0.001:
                    lamstr = r'\lambda/6'
                tline(d, '$Z_0=' + str(impedance) + unit + r',\, ' + lamstr + '$')
            if typus == "Open Stub":
                lamlen = ele.attr['length']
                lamstr = str(round(lamlen,3)) + r'\lambda'
                if abs(lamlen - 1. / 8) < 0.001:
                    lamstr = r'\lambda/8'
                if abs(lamlen - 1. / 4) < 0.001:
                    lamstr = r'\lambda/4'
                if abs(lamlen - 1. / 3) < 0.001:
                    lamstr = r'\lambda/3'
                if abs(lamlen - 1. / 6) < 0.001:
                    lamstr = r'\lambda/6'
                tstub(d, '$' + str(ele.value) + unit + '$', '$' + lamstr + '$')
            if self.showArrows:
                intermediateImpedanceLabel(d, "$Z_"+str(ele.id)+"$")
        return d



########################################################################################################################
class Smithpaper(Smith):

    def __init__(self, ax, typ='smith', Z0=1, fontscale=1, alpha= 1, fineness=1, showrings=False, **kwargs):
        self.f = 0
        self.ax = ax
        self.Z0 = Z0
        self.kwargs = kwargs
        self.linefactor = 1
        self.elements = []
        self.impedances = []
        self.elementvalues = []
        self.elementvaluesunit = []
        self.circuittypes = []
        self.circuit = []
        self.Zinlist = []
        self.isnorm = False
        self.paper = False
        self.showArrows = False
        self.linewidth = 1.0
        self.fontscale = fontscale
        self.resolution = 300
        self.fatfactor = 2.2
        self.show_neg_numbers = True
        e.style(e.STYLE_IEC)
        self.alpha = alpha
        
        if typ =="smith": 
            self.addimpedancegrid(color="#FF8888", alpha = alpha)
        elif typ =="inverted":
            self.addadmittancegrid(color="#8888FF", alpha = alpha)
        elif typ == "both":
            self.addadmittancegrid(color="#8888FF", alpha = alpha)
            self.addimpedancegrid(color="#FF8888", alpha = alpha)
            
            
        if showrings:
            self.addGammaAngleRing()
            self.addLambdaRing()
            self.addLambdaRing2()
            self.addOuterRing()


    def addimpedancegrid(self, inverted=False, fine = False, **kwargs):
        '''
        :param inverted: plots an inverted chart, aka Admittance Chart
        :param kwargs: plot arguments, e.g. color, lw, ...
        :return:
        '''
        color = "gray"
        if self.show_neg_numbers:
            labfront = "-"
        else:
            labfront = ""

        ## Cirlces ###########
        clist = []
        ## a group of circles consist start imag, end imag, resolution and Fat line yes/no
        an = [0,     0.2,   0.5,   1,    2,    5,   10,   20,   50]
        res = [ 0.01,  0.02,  .05,   .1,   0.2,   1,    2,   10  ]
        minor_ranges = [ [an[i], an[i+1] ] for i in range(len(an)-1) ]
        minor_res    = res
        extendto = 3
        labelmarker = (list(around(arange(0.1,1,0.1),1)))
        labelmarker += (list(around(arange(1, 2, 0.2), 1)))
        labelmarker += range(2,6)
        labelmarker += (10,20,50)
        llw = [0.4 * self.linefactor, 0.4 * self.linefactor * self.fatfactor]
        ## Outer Circle
        e1 = patches.Circle((0, 0), 1,linewidth=llw[1], fill=False)
        self.ax.add_patch(e1)

        ################################################################################################################
        ### Constant real part circles #################################################################################
        ################################################################################################################
        print("------------------")
        maxbox = len(minor_ranges)
        for seg1 in range(maxbox):
          for seg2 in range(maxbox):
            istart,istop  = minor_ranges[seg2]
            if seg1 == seg2: reso =  minor_res[seg2]
            elif seg1 > seg2: reso =  minor_res[seg1]
            elif seg1 < seg2: reso =  minor_res[seg2]
            for zreal in arange(*minor_ranges[seg1], reso):

                ## Find Center and radius of circles/segments
                ## see https://www.allaboutcircuits.com/technical-articles/mathematical-construction-and-properties-of-the-smith-chart/
                radius = 1 / (zreal + 1)
                center = zreal / (zreal + 1)
                Zstart = zreal + istart * 1j
                gammastart = (Zstart - 1) / (Zstart + 1)
                #ax.plot(real(gammastart), imag(gammastart) , "go")
                angle1 = arctan2(imag(gammastart) / radius, real(gammastart- center) / radius) * 180 / pi
                Zend = zreal + istop * 1j
                gammaend = (Zend - 1) / (Zend + 1)
                #ax.plot(real(gammaend), imag(gammaend), "ro")
                angle2 = arctan2(imag(gammaend) / radius, real(gammaend- center) / radius) * 180 / pi
                if inverted:
                    center = -center
                    angle1, angle2 = angle1 +180, angle2 +180
                e1 = patches.Arc((center, 0), 2 * radius, 2 * radius, linewidth=llw[0], theta1 = angle2,
                                 theta2 = angle1, **kwargs)
                e2 = patches.Arc((center, 0), 2 * radius, 2 * radius, linewidth=llw[0], theta1=-angle1,
                                 theta2=-angle2, **kwargs)
                self.ax.add_patch(e1)
                self.ax.add_patch(e2)

        ### Labels ##################

        for ll in labelmarker:
        #for ll in []:
            zreal = ll
            gamlabel = (zreal - 1.0) / (zreal + 1.0)
            orient = 1
            labeloffsets = [0.006,1.01,0.99]
            if inverted:
                 gamlabel = -gamlabel
                 orient = -1
                 labeloffsets = [-0.006,0.99,1.013]

            self.ax.text(gamlabel, labeloffsets[0], str(ll), rotation_mode='anchor', ha='left', va='bottom',
                          rotation=90 * orient, fontsize=6 * self.fontscale, **kwargs)
            ## Fat Line
            cend = an[extendto]
            if ll >= 0.5:  cend = an[extendto + 1]
            if ll > 2:  cend = an[extendto + 2]
            radius = 1 / (zreal + 1)
            center = zreal / (zreal + 1)
            Zend = zreal + cend * 1j
            gammaend = (Zend - 1) / (Zend + 1)
            angle2 = arctan2(imag(gammaend) / radius, real(gammaend - center) / radius) * 180 / pi
            angle1 = 180
            if inverted:
                center = -center
                angle1, angle2 = angle1 +180, angle2 +180
            e1 = patches.Arc((center, 0), 2 * radius, 2 * radius, linewidth=llw[1], theta1=angle2,theta2=angle1, **kwargs)
            e2 = patches.Arc((center, 0), 2 * radius, 2 * radius, linewidth=llw[1], theta1=angle1,theta2=-angle2, **kwargs)
            self.ax.add_patch(e1)
            self.ax.add_patch(e2)

            ### Rotated Labels ###
            if ll<= 1 and not round(ll*10) % 2:
                gamlabel = (zreal + labeloffsets[1]*1j - 1) / (zreal + labeloffsets[1]*1j + 1)
                anglabel = 180 + angle(gamlabel - 1 - 1j) * 180 / pi
                if inverted: gamlabel = -gamlabel
                self.ax.text(real(gamlabel), imag(gamlabel), str(ll), rotation_mode='anchor', ha='left',
                              va='bottom',
                              rotation=anglabel, fontsize=6 * self.fontscale,  **kwargs)
                gamlabel = (zreal + labeloffsets[2]*1j - 1) / (zreal + labeloffsets[2]*1j + 1)
                anglabel = 180 + angle(gamlabel - 1 - 1j) * 180 / pi
                if inverted: gamlabel = -gamlabel
                self.ax.text(real(gamlabel), imag(-gamlabel), str(ll), rotation_mode='anchor', ha='left',
                              va='bottom',
                              rotation=180 - anglabel, fontsize=6 * self.fontscale, **kwargs)
        self.ax.set_xlim(-1.2, 1.2)
        self.ax.set_ylim(-1.2, 1.2)
        self.ax.axis('off')
        self.ax.axis('equal')

        ################################################################################################################
        ## Constant imag Part Bends ####################################################################################
        ################################################################################################################

        ## horizontal real axis line ######################################
        self.ax.plot([-1,1], [0,0], linewidth=llw[1], **kwargs)
        #kwargs = {'color': "g"}
        print("------------------")
        for seg1 in range(maxbox):
            for seg2 in range(maxbox):
                istart, istop = minor_ranges[seg2]
                if seg1 == seg2:
                    reso = minor_res[seg2]
                elif seg1 > seg2:
                    reso = minor_res[seg1]
                elif seg1 < seg2:
                    reso = minor_res[seg2]
                myrange = minor_ranges[seg1]
                if (seg1,seg2) == (2,4):  myrange = [0.6,1] # unsolved yet
                if (seg1, seg2) == (1, 5):  myrange = [1, 1]  # unsolved yet
                if (seg1, seg2) == (2, 5):  myrange = [1, 1]  # unsolved yet
                if (seg1, seg2) == (1, 6):  myrange = [1, 1]  # unsolved yet
                for zimag in arange(*myrange, reso):
                    if zimag == 0: continue
                    radius = 1 / zimag
                    center = 1 / zimag
                    #ax.plot(1, center, "bo")
                    Zstart = istart + zimag *1j
                    gammastart = (Zstart - 1) / (Zstart + 1)
                    angle1 = arctan2(imag(gammastart-center*1j) / radius, real(gammastart-1) / radius) * 180 / pi
                    Zend = istop + zimag * 1j
                    gammaend = (Zend - 1) / (Zend + 1)
                    angle2 = arctan2(imag(gammaend-center*1j) / radius, real(gammaend-1) / radius) * 180 / pi
                    centerx = 1
                    if inverted:
                        center = -center
                        centerx = -centerx
                        angle1, angle2 = angle1 + 180, angle2 + 180
                    e1 = patches.Arc((centerx, center), 2 * radius, 2 * radius, linewidth=llw[0], theta1=angle1,
                                     theta2=angle2, **kwargs)
                    e2 = patches.Arc((centerx, -center), 2 * radius, 2 * radius, linewidth=llw[0], theta1=360-angle2,
                                     theta2=360-angle1, **kwargs)
                    self.ax.add_patch(e1)
                    self.ax.add_patch(e2)

        ### Labels imag part ##################
        for ll in labelmarker:
        #for ll in []:
            if ll == 0: continue
            zimag = ll
            gg = (1j * ll - 1) / (1j * ll + 1)
            if inverted: gg = -gg
            rad = 0.994
            x,y = rad * real(gg), rad * imag(gg)
            ang = angle(gg) * 180.0 / pi
            self.ax.text(x, y, str(ll), ha='right', va='bottom', rotation_mode='anchor', rotation=ang,
                              fontsize=6 * self.fontscale,  **kwargs)
            self.ax.text(x, -y, labfront + str(ll), ha='left', va='bottom', rotation_mode='anchor', rotation=180 - ang,
                              fontsize=6 * self.fontscale,  **kwargs)

            ## Fat Line imag bends
            cend = an[extendto]
            if ll > 1:  cend = an[extendto + 1]
            if ll > 2:  cend = an[extendto + 2]
            radius = 1 / zimag
            center = 1 / zimag

            Zstart = 0 + zimag * 1j
            gammastart = (Zstart - 1) / (Zstart + 1)
            #ax.plot(real(gammastart), imag(gammastart) , "go")
            angle1 = arctan2(imag(gammastart - center * 1j) / radius, real(gammastart - 1) / radius) * 180 / pi
            Zend = cend + zimag * 1j
            gammaend = (Zend - 1) / (Zend + 1)
            #ax.plot(real(gammaend), imag(gammaend) ,"ro")
            angle2 = arctan2(imag(gammaend - center * 1j) / radius, real(gammaend - 1) / radius) * 180 / pi
            centerx = 1
            if inverted:
                center = -center
                centerx = -centerx
                angle1, angle2 = angle1 + 180, angle2 + 180
            e1 = patches.Arc((centerx, center), 2 * radius, 2 * radius, linewidth=llw[1], theta1=angle1,
                             theta2=angle2, **kwargs)
            e2 = patches.Arc((centerx, -center), 2 * radius, 2 * radius, linewidth=llw[1], theta1=360 - angle2,
                             theta2=360 - angle1, **kwargs)
            self.ax.add_patch(e1)
            self.ax.add_patch(e2)

    def addadmittancegrid(self, **kwargs):
        self.addimpedancegrid(inverted=True, **kwargs)


    def addanglering(self, radius_inner = 1.03, radius_outer = 1.05,
                     major = arange(0, 361, 10), majorlabels = arange(0, 361, 10), minor = arange(0, 361, 10),
                     angleoffset = 0, **kwargs):
        '''
        Adds a ring with angles
        :param kwargs:
        :return:
        '''
        if not "color" in kwargs: kwargs['color'] = 'k'
        if not "lw" in kwargs: kwargs['lw'] = 1
        if not "alpha" in kwargs: kwargs['alpha'] = 1.0
        #print(f' Kwargs: {kwargs}')


        e1 = patches.Circle((0,0),radius_inner, fill=False, **kwargs)
        self.ax.add_patch(e1)

        for i in range(len(major)):
            tic = major[i]
            lab = str(majorlabels[i])
            ang = (tic + angleoffset) * pi / 180
            x = cos(ang)
            y = sin(ang)
            self.ax.plot((radius_inner * x, radius_outer * x), (radius_inner * y, radius_outer * y), **kwargs)
            self.ax.text(x * radius_outer, y * radius_outer, lab, ha='center', va='bottom', rotation_mode='anchor',
                         rotation= -90 + ang * 180.0 / pi,
                         fontsize= 5 * self.fontscale, color=kwargs['color'], alpha=kwargs['alpha'])

        kwargs['lw'] = 0.5 * kwargs['lw']
        for tic in minor:
            ang = tic * pi / 180
            x = cos(ang)
            y = sin(ang)

            self.ax.plot((radius_inner * x, (radius_outer-0.01) * x), (radius_inner * y, (radius_outer-0.01) * y),  **kwargs)

        self.ax.axis('equal')
        plt.tight_layout(pad=.0)

    def addGammaAngleRing(self, resolution=10):
        """
        Adds a ring with the angles of the reflection coefficient
        :param resolution:
        :return:
        """
        major = concatenate((arange(0, 181, resolution), arange(0, -180, -resolution)))
        minor = arange(0, 361, 2)
        self.addanglering(radius_inner=1.03, radius_outer=1.05, major=major, majorlabels=major, minor=minor)

    def addLambdaRing(self, resolution=0.01):
        """
        Adds a ring with the angles of the reflection coefficient
        :param resolution:
        :return:
        """

        labels = [f"{x:4.2g}$\\lambda$" for x in arange(0, 0.5, resolution)]
        major = arange(180, -180, -resolution*720)
        minor = arange(0, 360, resolution*720/5)
        self.addanglering(radius_inner=1.08, radius_outer=1.10, major=major, majorlabels= labels , minor=minor)

    def addLambdaRing2(self, resolution=0.01):
        """
        Adds a ring with the angles of the reflection coefficient
        towards generator
        :param resolution:
        :return:
        """

        labels = [f"{x:4.2g}$\\lambda$" for x in arange(0, 0.5, resolution)]
        major = arange(-180, 180, resolution*720)
        minor = arange(0, 360, resolution*720/5)
        self.addanglering(radius_inner=1.13, radius_outer=1.15, major=major, majorlabels= labels , minor=minor)

    def addOuterRing(self):
        self.addanglering(radius_inner=1.18, radius_outer=1.20, major=[], majorlabels=[], minor=[])


#######################################################################################
###Draw Schematics ####################################################################
#######################################################################################
#
# This plots the schematic of the circuit
# it uses a very old syntax of schemdraw, not recommended for new designs !!
#
#
from schemdraw.segments import *
    
class Transline(e.Element2Term):

    def __init__(self, l=3, *args, **kwargs):
        super().__init__(*args, **kwargs)
        lw = 2.8 * l/4
        lh = 0.25
        self.segments.append(Segment([[0, 0],[0,lh/2],[lw,lh/2],[lw,-lh/2], [0,-lh/2],[0,0],[lw,0]], fill= 'black'))


class Tstub(e.Element2Term):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        lw = 2.8
        lh = 0.25
        perspec = 2.0
        self.segments.append(
            Segment([[0, 0], [0, lh / 2], [lw, lh / 2 + perspec], [lw, -lh / 2 + perspec], [0, -lh / 2], [0, 0], [lw, 0]], fill='black'))


class ImpArrow(e.Element):

    def __init__(self, *d, ofst=(0,0), **kwargs):
        super().__init__(*d, **kwargs)
        lw = 0.5
        lh = 1.5
        a = 0.15
        off = 2.3
        x = ofst[0]
        y = ofst[1]
        self.segments.append(Segment([[0+x, off+y], [lw+x, off+y], [lw+x, off+lh+y]]) )
        self.segments.append(Segment([[0+x, off+y], [a+x, off+0.4*a+y], [a+x, off-0.4*a+y],[0+x,off+y]], fill='green') ) # arrowhead


def startR(d, lab='', lsh = 0.5):
    d.push()
    d += e.Line(l=lsh).right()
    d += e.Resistor().down().label(lab,'bot')
    d += e.Line(l=lsh).left()
    d.pop()


def startL(d, lab=''):
    d.push()
    d += e.Inductor().down().label(lab,'bot')
    d.pop()


def startC(d, lab=''):
    d.push()
    d += e.Capacitor.down().label(lab,'bot')
    d.pop()

def startport(d, lab="$P2$"):
    d.push()
    d += e.Dot(open=True)
    d += e.Gap().down().label(lab,'bot')
    d += e.Dot(open=True)
    d.pop()
    

def seriesR(d, lab='', lsh = 2):
    d.push()
    d += e.Resistor(l=lsh, flip=True).left().label(lab)
    d.pop()
    d += e.Gap().label('').down()
    d += e.Line(l=lsh).left()
    d += e.Gap().label('').up()


def seriesC(d, lab='', lsh = 2):
    d.push()
    d += e.Capacitor(l=lsh).left().label(lab)
    d.pop()
    d += e.Gap().label('').down()
    d += e.Line(l=lsh).left()
    d += e.Gap().label('').up()


def seriesL(d, lab='', lsh = 2):
    d.push()
    d += e.Inductor(l=lsh, flip=True).left().label(lab)
    d.pop()
    d += e.Gap().label('').down()
    d += e.Line(l=lsh).left()
    d += e.Gap().label('').up()

def seriesSpace(d, lab='', lsh = 2):
    d.push()
    d += e.Line(l=lsh).left().label(lab)
    d.pop()
    d += e.Gap().label('').down()
    d += e.Line(l=lsh).left()
    d += e.Gap().label('').up()

def shuntR(d, lab='', lsh=1):
    d += e.Line(l=lsh).left()
    d += e.Dot()
    d.push()
    d += e.Resistor().down().label(lab)
    d += e.Dot()
    d.push()
    d += e.Line(l=lsh).right()
    d.pop()
    d += e.Line(l=lsh).left()
    # d += e.DOT)
    d.pop()

    d += e.Line(l=lsh).left()
    # d += e.DOT)


def shuntC(d, lab='', lsh=1):
    d += e.Line(l=lsh).left()
    d += e.Dot()
    d.push()
    d += e.Capacitor().down().label(lab)
    d += e.Dot()
    d.push()
    d += e.Line(l=lsh).right()
    d.pop()
    d += e.Line(l=lsh).left()
    # d += e.DOT)
    d.pop()
    d += e.Line(l=lsh).left()
    # d += e.DOT)


def shuntL(d, lab='', lsh = 1):
    d += e.Line(l=1.).left()
    d += e.Dot()
    d.push()
    d += e.Inductor().down().label(lab)
    d += e.Dot()
    d.push()
    d += e.Line(l=lsh).right()
    d.pop()
    d += e.Line(l=lsh).left()
    # d += e.DOT)
    d.pop()
    d += e.Line(l=lsh).left()
    # d += e.DOT)


def tline(d, lab="$Z_0$", l=3):
    d += e.Line(l=0.5).left()
    if hide: lab = "$Z_0$"
    d += Transline(l=l).left().label(lab)
    d += e.Line(l=0.5).left()
    d.push()
    d += e.Gap().label('').down()
    d += e.Line(l=0.5).right()
    d += Transline(l=l).right()
    d += e.Line(l=0.5).right()
    d.pop()
    # d += e.DOT)


def tstub(d, lab="$Z_0$", lab2=''):
    d += e.Line(l=0.9).left()
    d.push()
    d.push()
    if hide: lab = "     $Z_0$"; lab2 = ""
    d += Tstub().right().label(lab+' '+lab2, ofst=(-1.1,-1.6),rotate=35)
    d.pop()
    d += e.Gap().label('').down()
    d.push()
    d += e.Line(l=0.9).right()
    d.pop()
    d.push()
    d += Tstub.right()
    d.pop()
    d += e.Line(l=0.5).left()
    d.pop()

def inputport(d, lab=""):
    d.push()
    d += e.Line(l=0.9).left()
    d += e.Dot(open=True)
    d += e.Gap().down().label(lab,'top')
    d += e.Dot(open=True)
    d += e.Line(l=0.9).right()
    d.pop()
    
def generator(d, lab="$U_0$",labR="$R_s$"):
    d.push()
    d += e.Dot(open=True)
    d += e.Line(l=0.3).left()
    d += e.Resistor().left().label(labR)
    d += e.SourceSin().down().label(lab)
    d += e.Line(l=3.3).right()
    d += e.Dot(open=True)
    d.pop()


def portannotation(d,lab,**kwargs):
    d.push()
    d += e.Line(l=0.45).left()
    d += (mygap:= e.Gap().down().label(""))
    d += e.CurrentLabel().at(mygap).down().label(lab)
    d += e.Line(l=0.45).right()
    d.pop()

def intermediateImpedanceLabel(d,lab="$Z$",**kwargs):
    d.push()
    d += ImpArrow(color="green").label(label=lab,loc="bottom",ofst=(0.3,0.2)).left()
    d.pop()


def twoport(d, lab = "[S]"):
    TwoPort = e.Ic(pins=[e.IcPin(side='left'),
                  e.IcPin(side='left', anchorname='Port1'),
                  e.IcPin(side='right'),
                  e.IcPin(side='right', anchorname='Port2')],
            edgepadW = 1.0,  # Make it a bit wider
            edgepadH = 1.0,
            pinspacing=3).label(lab, 'center', fontsize=28)
    x,y = d.here   
    d += TwoPort.right().anchor('Port2') 
    d.move(x-2,y)   

class smith(Smith):  # Definition for historical reasons when using small smith class def
    pass

def lSectionMatch(Zl,Z0=50):
    '''
    Solve the matching network as lumped L-section following Pozar page 224
    :param Zl:
    :param Z0:
    :return: The L-Section reactance X and susceptance B, (Series element left? boolean)
    '''
    R = real(Zl)
    X = imag(Zl)

    if R>Z0:
        B1 = ( X + sqrt(R/Z0) * sqrt(R**2+X**2-Z0*R) ) / (R**2+X**2)
        X1 = 1/B1 + X*Z0/R - Z0/B1/R
        B2 = (X - sqrt(R / Z0) * sqrt(R ** 2 + X ** 2 - Z0 * R)) / (R ** 2 + X ** 2)
        X2 = 1 / B2 + X * Z0 / R - Z0 / B2 / R
        serleft = True
    else:
        X1 = sqrt(R*(Z0-R))-X
        B1 = sqrt((Z0-R)/R) / Z0
        X2 = -sqrt(R * (Z0 - R)) - X
        B2 = -sqrt((Z0 - R) / R) / Z0
        serleft = False

    return (X1,B1),(X2,B2), serleft

def stubLineMatch(Zl,Z0=50):
    '''
    Solve the matching network as line Network with open Stub
    :param Zl: load impedance (complex)
    :param Z0: line and input impedance
    :return: l1,l2 (The length of the stub l2 and the inserted line l1)
    '''

    Gammamatch = (Zl-Z0) / (Zl+Z0)
    Gammamatch = conj(Gammamatch) # since we start at load
    Z1 = Z0

    ### finding length of open stub line
    for len2 in arange(0.0001, 1, 0.0001):
        betal = 2 * pi * len2
        Zstub = -1j * Z0 / tan(betal)
        Z2 = Z1 * Zstub / (Z1 + Zstub)
        Gam = (Z2 - Z0) / (Z2 + Z0)
        if abs(Gam) > abs(Gammamatch): break
    ### finding length of inserted line
    if (angle(Gam) > angle(Gammamatch)):
        betal = (angle(Gam) - angle(Gammamatch)) / 2
    else:
        betal = (angle(Gam) - angle(Gammamatch)) / 2 + pi
    len1 = round(betal / 2 / pi, 3)
    return len1, len2

########################################################################################################################



if __name__ == "__main__":
    # import doctest
    # doctest.testmod()

    demo = 5

    if demo == -1:
        ## Plot Empty smithchart
        fig, ax = plt.subplots(figsize=(16,16))
        Z0 = 50
        mysmith = Smith(ax, 'both', fineness= 3)
        mysmith.addpolargrid(alpha= 0.3)
        mysmith.addanglering(alpha=0.6, color = "g")
        #mysmith = Smith(ax, 'smith', Z0, color = "g", lw = 1, alpha = 0.5)
        plt.show()

    if demo == 0:
        fig, ax = plt.subplots()
        Z0 = 50
        mysmith = Smith(ax, 'both', Z0)
        Z1 = mysmith.addstart(20 )
        Z2 = mysmith.addline(Z1, 0.1)
        Z3 = mysmith.addpara(Z2,-20j)
        mysmith.addangle(Z2)
        mysmith.addangle(Z3)
        Z4 = mysmith.addseries(Z3, -60j)
        mysmith.addpoint(Z1, '$Z_1$', 'SW')
        mysmith.addpoint(Z2, '$Z_2$', 'NE')
        mysmith.addpoint(Z3, '$Z_3$', 'NE')
        mysmith.addpoint(Z4, '$Z_in$', 'SE')
        plt.show()

    if demo == 5:
        ## Plot smithchart in real Smith Chart A4
        fig, ax = plt.subplots(figsize=(8,8))
        Z0 = 50
        mysmith = Smithpaper(ax, 'smith', Z0, showrings=True)

        if True:
            Z1 = mysmith.addstart(20 )
            Z2 = mysmith.addline(Z1, 0.1)
            Z3 = mysmith.addpara(Z2,-20j)
            mysmith.addangle(Z2)
            mysmith.addangle(Z3)
            Z4 = mysmith.addseries(Z3, -60j)
            mysmith.addpoint(Z1, '$Z_1$', 'SW')
            mysmith.addpoint(Z2, '$Z_2$', 'NE')
            mysmith.addpoint(Z3, '$Z_3$', 'NE')
            mysmith.addpoint(Z4, '$Z_{in}$', 'SE')
            #mysmith.addrulermarker(0.7)
        plt.tight_layout()
        plt.savefig("smith.svg")
        plt.savefig("smith.png",dpi=300)
        plt.show()

    #### Complex Demo with circuit schematic generation
    if demo == 1:
        ## Plot smithchart
        fig, ax = plt.subplots(figsize=(10,10))
        Z0 = 50
        mysmith = Smith(ax, 'both', Z0, fineness=3)
        # mysmith.addpolargrid()
        mysmith.showArrows = True
        Z1 = mysmith.addstart(20 - 10j)
        Z2 = mysmith.addpara(Z1, 30j)
        Z3 = mysmith.addline(Z2, 0.3)
        mysmith.addangle(Z2)
        mysmith.addangle(Z3)
        Z4 = mysmith.addseries(Z3, -60j)
        mysmith.addinput()
        mysmith.addpoint(Z1, '$Z_1$', 'SW')
        mysmith.addpoint(Z2, '$Z_2$', 'NE')
        mysmith.addpoint(Z3, '$Z_3$', 'SW')
        mysmith.addpoint(Z4, '$Z_in$', 'SE')

        for ele in mysmith.circuit:
            Z = ele.Zin
            n = ele.id
            Gam = magphasetex((Z - Z0) / (Z + Z0))
            print('$Z_{0:1} = ({1:4.2f}) \\;\\Omega \\qquad \\Gamma_{0:1} = {2:s}$ \n'.format(n, Z, Gam))

        Zin = Z4
        gam = (Zin - Z0) / (Zin + Z0)
        mysmith.addarrow(gam, color='b', lw=2)
        mysmith.addruler()
        mysmith.addrulerdistance(gam)
        mysmith.addrulermarker(gam, color="green")
        plt.show()

        ### Create a circuit schematic of the analysed circuit  ##
        d = mysmith.plotschematic()
        d.draw()
        plt.show()

    if demo == 2:
        ### Example of circles #####################################
        fig, ax = plt.subplots(figsize=(10, 10))
        Z0 = 50
        cols = ['k', 'red', 'green', 'yellow', 'powderblue']
        mysmith = Smith(ax, 'both', Z0, fineness=2)
        for i in range(1, 5):
            C = (0.2 + 0.5j) + i * 0.1 * (0.3 + 0.3j)
            R = 0.4 - i * 0.1
            mysmith.addcircle(C, R, color=cols[i])
        plt.show()
        plt.close(fig)
        print("Done")

    ### Frequency Sweep of C-LCR Circuit
    if demo == 3:
        ### Example of circles #####################################
        fig, ax = plt.subplots(figsize=(10, 10))
        Z0 = 50
        f = arange(.1e9, 5e9, .001e9)
        omega = 2 * pi * f
        C = 10e-12
        L = 10e-9
        R = 30
        C2 = 40e-12
        Yc = 1j * omega * C
        Yl = 1 / (1j * omega * L)
        Yr = 1 / R
        Ypar = Yc + Yl + Yr
        Zpar = 1 / Ypar
        Z = Zpar + 1 / (1j * omega * C2)
        Gam = (Z - Z0) / (Z + Z0)

        mysmith = Smith(ax, 'paper', Z0, fineness=1)
        ax.plot(real(Gam), imag(Gam), lw=2)
        #ax.set_title('Frequency Sweep of C-LRC Circuit', fontsize=20)

        ### find resonance and annotate ##############################
        omega0 = sqrt(1 / L / C)
        fres = omega0 / 2 / pi
        idx = argmin(abs(f - fres))
        x, y = (real(Gam[idx]), imag(Gam[idx]))
        ax.plot(x, y, 'ko')
        ftext = "$f_{{res}} ={0:5.2f} $ MHz".format(fres / 1e6)
        arrow_properties = dict(facecolor="black", width=0.5, headwidth=4, shrink=0.04)
        plt.annotate(
            ftext, xy=(x, y), xytext=(0.1, 0.2), arrowprops=arrow_properties, fontsize=16)
        plt.savefig("out.svg")
        #plt.show()
        file = superimposeChart(fig)
        plt.close(fig)
        print("Done")

        ### L syle Line Matching
        if demo == 4:
            Z0 = 100
            Zl = 25 - 200j  # invent a point !
            sol1, sol2, serleft = lSectionMatch(Zl, Z0)
            print(sol1, sol2)
            for sol in (sol1, sol2):
                X, B = sol
                fig, ax = plt.subplots(figsize=(8, 8))
                mysmith = Smith(ax, 'both', Z0)
                Z1 = mysmith.addstart(Zl)
                if serleft:
                    Z2 = mysmith.addpara(Z1, 1 / B / 1j)
                    Z3 = mysmith.addseries(Z2, X * 1j)
                else:
                    Z2 = mysmith.addseries(Z1, X * 1j)
                    Z3 = mysmith.addpara(Z2, 1 / B / 1j)
                mysmith.addpoint(Z1, '$Z_1$', 'SW')
                mysmith.addpoint(Z2, '$Z_2$', 'NE')
                mysmith.addpoint(Z3, '$Z_{in}', 'SE')
                # plt.show()

            l1, l2 = stubLineMatch(Zl, Z0)
            print(l1, l2)
            fig, ax = plt.subplots(figsize=(8, 8))
            mysmith = Smith(ax, 'both', Z0)
            mysmith.addanglering()
            Z1 = mysmith.addstart(Zl)
            Z2 = mysmith.addline(Z1, l1, Z0)
            Z3 = mysmith.addstubopen(Z2, l2, Z0)
            mysmith.addangle(Z1)
            mysmith.addangle(Z2)
            mysmith.addpoint(Z1, '$Z_1$', 'SW')
            mysmith.addpoint(Z2, '$Z_2$', 'NE')
            mysmith.addpoint(Z3, '$Z_{in}$', 'SE')

            Zopen = 1e55
            Zopen = mysmith.addstart(Zopen)
            Zstub = mysmith.addline(Zopen, l2, Z0)
            mysmith.addangle(Zstub)
            mysmith.addangle(Zopen)
            mysmith.addpoint(Zopen, 'Open', 'SE')
            mysmith.addpoint(Zstub, '$Z_{stub}$', 'SE')

            plt.show()

        #### Cartoon Style ##################################################
        if False:
            plt.xkcd()
            fig, ax = plt.subplots(figsize=(10, 10))
            mysmith = Smith(ax, 'smith', Z0, fineness=1)
            ax.plot(real(Gam), imag(Gam), lw=2)
            ax.set_title('Frequency Sweep of C-LRC Circuit', fontsize=20)
            ax.plot(x, y, 'ko')
            plt.annotate(
                ftext, xy=(x, y), xytext=(0.1, 0.2), arrowprops=arrow_properties, fontsize=16)
            plt.show()
            plt.close(fig)
            print("Done")
