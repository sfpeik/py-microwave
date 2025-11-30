"""
MoM (Modern Python) - Pocklington & Hallen implementations (variant B)

inspired from Balanis MoM Program 

Requirements:
    numpy, scipy, matplotlib

Usage examples at bottom.
"""
import numpy as np
import math
import cmath
from scipy.linalg import toeplitz, solve, lu_factor, lu_solve
import matplotlib.pyplot as plt
from typing import Tuple, Optional


# ---------- Utilities ----------
PI = math.pi
RAD = PI / 180.0
BETA = 2.0 * PI   # wave number in normalized (wavelength) units: beta = 2*pi / lambda; for lambda=1 => 2*pi
ETA = 120.0 * PI


def _ensure_odd(n: int) -> int:
    return n if n % 2 == 1 else n + 1


# ---------- Pocklington (pulse basis, point matching) ----------
def solve_pocklington(nm: int, tl: float, ra: float, iex: int = 2,
                      use_files: bool = True, plot_results: bool = True
                      ) -> Tuple[np.ndarray, complex, np.ndarray]:
    """
    Solve Pocklington formulation (pulse expansion, delta-gap or magnetic frill).
    Parameters:
        nm: number of subdivisions (should be odd)
        tl: total dipole length (in wavelengths)
        ra: wire radius (in wavelengths)
        iex: excitation type: 1=magnetic-frill, 2=delta-gap
        use_files: whether to write Curr-MoM_m.dat and Patt-MoM_m.dat
        plot_results: show matplotlib plots
    Returns:
        currents_half : complex array of currents on half dipole (length (nm+1)/2)
        zin : complex input impedance
        pattern_db : 1D array of pattern (degrees 0..180) in dB
    """
    nm = _ensure_odd(nm)
    hl = tl / 2.0
    nmh = int(0.5 * (nm + 1))
    dz = 2.0 * hl / nm
    # Simpson integration parameters
    n = 79
    b = 0.5 * dz
    a = -0.5 * dz
    hd = (b - a) / (n + 1)

    # Preallocate vector a (toeplitz first row+col form in MATLAB version)
    # We'll compute kernel entries for center-to-center interaction (Toeplitz property)
    zmn = np.zeros(2 * nm - 1, dtype=complex)

    for I in range(1, nm + 1):
        zn = hl - (I - 0.5) * dz
        za1 = zn - (hl - 0.5 * dz) + a  # but original code used zm = hl - 0.5*dz: careful -> easier computing directly
        # We'll follow MATLAB more directly: compute integral via Simpson on interval [a,b] for each segment center difference
        # Compute Simpson endpoints:
        # cgp1 at x=a, cgp2 at x=b, plus 79 internal nodes
        crt = 0+0j
        # compute for x = a, b and internal nodes
        for x, weight in [(a, 1.0), (b, 1.0)]:
            zx1 = zn - (hl - 0.5 * dz) + x  # but simpler: integrate variable x over [a,b] around center zm = hl - 0.5*dz
            # To faithfully reproduce original, compute r with radius ra
            # However original uses zm = hl - 0.5*dz constant, so use that:
            zm = hl - 0.5 * dz
            zx = zn - zm + x
            r = math.sqrt(ra * ra + zx * zx)
            # kernel expression from MATLAB
            val = cmath.exp(-1j * BETA * r) * ((1.0 + 1j * BETA * r) * (2.0 * r * r - 3.0 * ra * ra) + (BETA * ra * r) ** 2) / (2.0 * BETA * r ** 5)
            crt += weight * val

        # internal nodes
        for k in range(1, n + 1):
            xk = a + k * hd
            zx = zn - (hl - 0.5 * dz) + xk
            r = math.sqrt(ra * ra + zx * zx)
            val = cmath.exp(-1j * BETA * r) * ((1.0 + 1j * BETA * r) * (2.0 * r * r - 3.0 * ra * ra) + (BETA * ra * r) ** 2) / (2.0 * BETA * r ** 5)
            if (k % 2) != 0:
                crt += 4.0 * val
            else:
                crt += 2.0 * val

        crt = crt * hd * (1.0 / 3.0)
        zmn[I - 1] = crt
        if I != 1:
            zmn[nm + I - 2] = crt  # MATLAB placed at zmn(nm+I-1)

    # excitation vector cga
    cga = np.zeros(nm, dtype=complex)
    rb = 2.3 * ra
    tlab = 2.0 * math.log(2.3)
    for i in range(1, nm + 1):
        zi = hl - (i - 0.5) * dz
        r1 = BETA * math.sqrt(zi * zi + ra * ra)
        r2 = BETA * math.sqrt(zi * zi + rb * rb)
        if iex == 1:
            cga[i - 1] = -1j * (BETA ** 2) / (ETA * tlab) * (cmath.exp(-1j * r1) / r1 - cmath.exp(-1j * r2) / r2)
        else:
            if i != nmh:
                cga[i - 1] = 0.0 + 0.0j
            else:
                cga[i - 1] = -1j * BETA / (ETA * dz)

    # Build full Toeplitz matrix from zmn vector: MATLAB's zmn holds first row then first column.
    # First element a11 = zmn(1). First row = zmn[0:nm]. First column = [zmn[0], zmn[nm:]]? Let's reconstruct simpler:
    # We assume symmetric Toeplitz with entries t_k = zmn[k-1] for k=1..nm (distance multiple).
    # So build first column as t0, t1, ..., t_{nm-1}
    t_col = np.zeros(nm, dtype=complex)
    for k in range(nm):
        t_col[k] = zmn[k]  # matches MATLAB where first nm entries are first-row
    # form Toeplitz matrix
    A = toeplitz(t_col)
    # Solve A x = cga
    x = solve(A, cga)

    # currents for half dipole
    currents_half = x[:nmh]

    # write currents to file if requested
    if use_files:
        with open("Curr-MoM_m.dat", "w") as fid:
            fid.write("CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE\n")
            fid.write("POSITION Z:   MAGNITUDE:       REAL PART:       IMAGINARY PART:      SECTION #:\n\n")
            for i in range(1, nmh + 1):
                xi = hl - (i - 0.5) * dz
                yi = abs(currents_half[i - 1])
                fid.write(f"{xi:1.4f}       {yi:1.6f}            {currents_half[i - 1].real:1.6f}        {currents_half[i - 1].imag:1.6f}          {i:2d}\n\n")

    # plot current magnitude along half dipole
    if plot_results:
        i_idx = np.arange(1, nmh + 1)
        xi = hl - (i_idx - 0.5) * dz
        yi = np.abs(currents_half)
        plt.stem(xi, yi, basefmt=" ")
        plt.grid(True)
        plt.xlabel("Position (z)")
        plt.ylabel("Magnitude")
        plt.title("CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE (Pocklington)")
        plt.show()

    # input impedance
    zin = 1.0 / x[nmh - 1]

    # radiation pattern (0..180 deg)
    etmm = np.zeros(181)
    for i in range(1, 182):
        theta = (i - 1) * RAD
        cth = math.cos(theta)
        sth = math.sin(theta)
        if abs(cth) < 0.001:
            ft = 1.0
        else:
            ft = math.sin(BETA * dz * cth * 0.5) / (BETA * dz * cth * 0.5)
        crt = 0+0j
        for m in range(1, nm + 1):
            zm_loc = hl - (m - 0.5) * dz
            crt += cmath.exp(1j * BETA * zm_loc * cth) * ft * x[m - 1] * dz
        ptt = abs(crt) * (sth ** 2) * ETA * 0.5
        etmm[i - 1] = ptt

    # normalize and convert to dB
    amax = np.max(etmm)
    pattern_db = np.zeros_like(etmm)
    for i in range(len(etmm)):
        ptt = etmm[i] / amax
        if ptt <= 1e-5:
            ptt = 1e-5
        pattern_db[i] = 20.0 * math.log10(ptt)

    if use_files:
        with open("Patt-MoM_m.dat", "w") as pok:
            pok.write("RADIATION PATTERN   vs   OBSERVATION ANGLE THETA\n")
            pok.write("THETA (in degrees)        MAGNITUDE(in dB)\n")
            for i in range(1, 182):
                xi_deg = i - 1
                pok.write(f"{xi_deg:3.1f}                            {pattern_db[i - 1]:7.3f}\n")
            # mirror for 181..360 like original
            for i in range(182, 362):
                xi_deg = i - 1
                pok.write(f"{xi_deg:3.1f}                            {pattern_db[362 - i]:7.3f}\n")

    # polar plot (0..360)
    if plot_results:
        etmm1 = np.concatenate([pattern_db, pattern_db[-2::-1]])  # [0..180] + [179..0] -> 361 values
        angles = np.deg2rad(np.arange(0, 361))
        # convert dB to linear for polar plot radius
        linear = 10.0 ** (etmm1 / 20.0)
        ax = plt.subplot(projection='polar')
        ax.plot(angles, linear)
        ax.set_title("RADIATION PATTERN (Pocklington) - normalized")
        plt.show()

    return currents_half, zin, pattern_db


# ---------- Hallen (full matrix + Gaussian quadrature) ----------
def solve_hallen(nm: int, tl: float, ra: float,
                 use_files: bool = True, plot_results: bool = True) -> Tuple[np.ndarray, complex, np.ndarray]:
    """
    Solve Hallen's integral equation via matrix assembly and LU solve.
    Parameters:
        nm: number of subdivisions (odd or even)
        tl: total dipole length (wavelengths)
        ra: wire radius (wavelengths)
    Returns:
        elecur: complex current vector
        zin: complex input impedance
        pattern_db: radiation pattern (0..180 deg) in dB
    """
    n = nm
    l = tl
    rho = ra
    rtod = 180.0 / PI
    eta = ETA
    bk = BETA
    cj = 1j
    dz = l / (2.0 * (n - 1))   # matches MATLAB
    dz1 = l / (2.0 * n)

    # 16-point Gaussian quadrature nodes & weights from MATLAB code:
    absica = np.array([
        -0.095012509837637, -0.281603550779259,
        -0.458016777657227, -0.617876244402644,
        -0.755404408355003, -0.865631202387832,
        -0.944575023073233, -0.989400934991650,
         0.095012509837637,  0.281603550779259,
         0.458016777657227,  0.617876244402644,
         0.755404408355003,  0.865631202387832,
         0.944575023073233,  0.989400934991650
    ])
    wght = np.array([
        0.189450610455068, 0.182603415044924,
        0.169156519395002, 0.149595988816577,
        0.124628971255534, 0.095158511682493,
        0.062253523938648, 0.027152459411754,
        0.189450610455068, 0.182603415044924,
        0.169156519395002, 0.149595988816577,
        0.124628971255534, 0.095158511682493,
        0.062253523938648, 0.027152459411754
    ])

    # Preallocate matrices
    zmatrx = np.zeros((n, n), dtype=complex)
    elecur = np.zeros(n, dtype=complex)

    # Fill matrix & excitation vector
    for i in range(1, n + 1):
        z = (2 * i - 1) * dz1 / 2.0
        zmatrx[i - 1, n - 1] = -math.cos(bk * z)  # MATLAB placed zmatrx(i,n) = -cos(bk*z)
        elecur[i - 1] = -cj * math.sin(bk * z) / (2.0 * eta)
        # now fill other columns
        for j in range(1, n):
            lower = (j - 1) * dz1
            upper = j * dz1
            no = 10
            delt = (upper - lower) / (2.0 * no)
            summation = 0+0j
            for t in range(1, no + 1):
                s = lower + (2 * t - 1) * delt
                for q in range(16):
                    x = s + absica[q] * delt
                    r1 = math.sqrt(rho * rho + (z - x) * (z - x))
                    r2 = math.sqrt(rho * rho + (z + x) * (z + x))
                    kernel = cmath.exp(-cj * bk * r1) / (4.0 * PI * r1) + cmath.exp(-cj * bk * r2) / (4.0 * PI * r2)
                    summation += wght[q] * kernel
            res = summation * delt
            zmatrx[i - 1, j - 1] = res

    # Scale pivoting info (like original)
    scal = np.zeros(n)
    for i in range(n):
        zmax = np.max(np.abs(zmatrx[i, :]))
        scal[i] = 1.0 / zmax if zmax != 0 else 1.0

    # LU with partial pivoting using scipy (safer)
    lu, piv = lu_factor(zmatrx)
    sol = lu_solve((lu, piv), elecur)

    elecur = sol
    # input impedance
    zin = 1.0 / elecur[0]

    # write Curr-MoM_m.dat
    if use_files:
        with open("Curr-MoM_m.dat", "w") as fid:
            fid.write("CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE\n")
            fid.write("POSITION Z:   CURRENT MAGNITUDE:    CURRENT PHASE:    SECTION:\n")
            for i in range(1, n):
                cur = abs(elecur[i - 1])
                pha = rtod * math.atan2(elecur[i - 1].imag, elecur[i - 1].real)
                fid.write(f"{(i * dz - dz / 2.0):2.6f}       {cur:2.6f}              {pha:3.6f}          {i}\n\n")

    # plot current magnitude
    if plot_results:
        idx = np.arange(1, n)
        cur = np.abs(elecur[:n - 1])
        plt.stem(idx * dz - dz / 2.0, cur, basefmt=" ")
        plt.grid(True)
        plt.xlabel("Position (z)")
        plt.ylabel("Current magnitude")
        plt.title("CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE (Hallen)")
        plt.show()

    # radiation pattern (0..pi)
    maxang = 181
    pwr = np.zeros(maxang)
    pmax = -1.0
    for j in range(1, maxang + 1):
        theta = (j - 1) / (maxang - 1) * PI
        sth = math.sin(theta)
        cth = math.cos(theta)
        arg = PI * dz * cth
        if abs(arg) <= 0.001:
            ft = 1.0
        else:
            ft = math.sin(arg) / arg
        crt = 0+0j
        n2 = n - 1
        for k in range(1, n2 + 1):
            argp = PI * (-l + (k - 1) * 2.0 * dz + dz) * cth
            crt += cmath.exp(1j * argp) * ft * elecur[k - 1]
            argp = PI * (-l + (n2 + k - 1) * 2.0 * dz + dz) * cth
            crt += cmath.exp(1j * argp) * ft * elecur[k - 1]
        power = abs(crt) * (sth ** 2)
        pwr[j - 1] = power
        if power >= pmax:
            pmax = power

    # convert to dB
    pattern_db = np.zeros(maxang)
    for i in range(maxang):
        pattrn = pwr[i] / pmax if pmax != 0 else 0.0
        if pattrn <= 1e-5:
            pattern_db[i] = -100.0
        else:
            pattern_db[i] = 20.0 * math.log10(pattrn)

    if use_files:
        with open("Patt-MoM_m.dat", "w") as hel:
            hel.write("RADIATION PATTERN   vs  OBSERVATION ANGLE THETA\n")
            hel.write("THETA (in degrees)      MAGNITUDE (in dB)\n")
            for i in range(1, maxang + 1):
                theta_deg = (i - 1) / (maxang - 1) * 180.0
                hel.write(f"{theta_deg:3.1f}                         {pattern_db[i - 1]:7.3f}\n")
            # mirror
            for i in range(182, 362):
                theta_deg = (i - 1) / (maxang - 1) * 180.0
                hel.write(f"{theta_deg:3.1f}                         {pattern_db[362 - i]:7.3f}\n")

    # polar plot
    if plot_results:
        patt = np.concatenate([pattern_db, pattern_db[-2::-1]])
        angles = np.deg2rad(np.arange(0, 361))
        linear = 10.0 ** (patt / 20.0)
        ax = plt.subplot(projection='polar')
        ax.plot(angles, linear)
        ax.set_title("RADIATION PATTERN (Hallen) - normalized")
        plt.show()

    return elecur, zin, pattern_db


# ---------- Usage examples ----------
if __name__ == "__main__":
    # Example: Pocklington
    nm_ex = 101  # odd
    tl_ex = 1.0  # wavelength units
    ra_ex = 0.001  # small radius in wavelengths
    print("Running Pocklington example...")
    currents_half, zin_p, patt_p = solve_pocklington(nm_ex, tl_ex, ra_ex, iex=2,
                                                     use_files=True, plot_results=True)
    print("Pocklington Zin =", zin_p)

    # Example: Hallen
    nm_h = 80
    tl_h = 1.0
    ra_h = 0.001
    print("Running Hallen example...")
    elecur, zin_h, patt_h = solve_hallen(nm_h, tl_h, ra_h, use_files=True, plot_results=True)
    print("Hallen Zin =", zin_h)

