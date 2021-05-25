#
# Vaious examples of the usage of the smith modules
#
#
import matplotlib.pyplot as plt
from numpy import *
from smith import *

demo = 1

if demo == 0:
    ## Plot smithchart
    fig, ax = plt.subplots(figsize=(16, 16))
    Z0 = 50
    mysmith = Smith(ax, 'both', Z0)
    Z1 = mysmith.addstart(20)
    Z2 = mysmith.addline(Z1, 0.1)
    Z3 = mysmith.addpara(Z2, -20j)
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
    fig, ax = plt.subplots()
    Z0 = 50
    mysmith = Smith(ax, 'paper', Z0)
    Z1 = mysmith.addstart(20)
    Z2 = mysmith.addline(Z1, 0.1)
    Z3 = mysmith.addpara(Z2, -20j)
    mysmith.addangle(Z2)
    mysmith.addangle(Z3)
    Z4 = mysmith.addseries(Z3, -60j)
    mysmith.addpoint(Z1, '$Z_1$', 'SW')
    mysmith.addpoint(Z2, '$Z_2$', 'NE')
    mysmith.addpoint(Z3, '$Z_3$', 'NE')
    mysmith.addpoint(Z4, '$Z_in$', 'SE')
    mysmith.addrulermarker(0.7)
    file = superimposeChart(fig)
    plt.close(fig)

#### Complex Demo with circuit schematic generation
if demo == 1:
    ## Plot smithchart
    fig, ax = plt.subplots(figsize=(10, 10))
    Z0 = 50
    mysmith = Smith(ax, 'both', Z0, fineness=1)
    # mysmith.addpolargrid()

    Z1 = mysmith.addstart(20 - 10j)
    Z2 = mysmith.addpara(Z1, 30j)
    Z3 = mysmith.addline(Z2, 0.3)
    mysmith.addangle(Z2)
    mysmith.addangle(Z3)
    Z4 = mysmith.addseries(Z3, -60j)
    mysmith.addpoint(Z1, '$Z_1$', 'SW')
    mysmith.addpoint(Z2, '$Z_2$', 'NE')
    mysmith.addpoint(Z3, '$Z_3$', 'SW')
    mysmith.addpoint(Z4, '$Z_in$', 'SE')

    for ele in mysmith.circuit:
        Z = ele.Zin
        n = ele.id
        Gam = magphasetex((Z - Z0) / (Z + Z0))
        print('$Z_{0:1} = ({1:4.2f}) \;\Omega \qquad \Gamma_{0:1} = {2:s}$ \n'.format(n, Z, Gam))

    Zin = Z4
    gam = (Zin - Z0) / (Zin + Z0)
    mysmith.addarrow(gam, color='b', lw=2)
    mysmith.addruler()
    mysmith.addrulerdistance(gam)
    mysmith.addrulermarker(gam, color="green")

    ### Create a circuit schematic of the analysed circuit  ##
    fig2, ax2 = plt.subplots(figsize=(7, 2))
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    d = mysmith.plotschematic()
    d.draw(ax2, show=False)  # must be set to false, otherwise plot will draw already in schemdraw module
    plt.tight_layout()
    plt.show()
    plt.close(fig)
    plt.close(fig2)

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
    # ax.set_title('Frequency Sweep of C-LRC Circuit', fontsize=20)

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
    # plt.show()
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
