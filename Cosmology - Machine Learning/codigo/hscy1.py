"""Utilidades para el catalogo de formas HSC-Y1 (S16A): carga, calibracion,
proyeccion a plano tangente, grilla de shear y Kaiser-Squires."""
import numpy as np
from astropy.io import fits

COLS = dict(w='ishape_hsm_regauss_derived_shape_weight',
            m='ishape_hsm_regauss_derived_shear_bias_m',
            c1='ishape_hsm_regauss_derived_shear_bias_c1',
            c2='ishape_hsm_regauss_derived_shear_bias_c2',
            sige='ishape_hsm_regauss_derived_sigma_e',
            erms='ishape_hsm_regauss_derived_rms_e',
            ra='ira', dec='idec',
            e1='ishape_hsm_regauss_e1', e2='ishape_hsm_regauss_e2',
            zp='photoz_best')
ZBINS = [0.3, 0.6, 0.9, 1.2, 1.5]

def cargar(path):
    d = fits.getdata(path, 1)
    ok = np.ones(len(d), bool)
    for c in COLS.values():
        ok &= ~d[c + '_isnull']
    cat = {k: np.asarray(d[c][ok], dtype=np.float64) for k, c in COLS.items()}
    ok2 = np.isfinite(cat['zp']) & np.isfinite(cat['e1']) & (cat['w'] > 0)
    return {k: v[ok2] for k, v in cat.items()}

def bin_tomografico(cat, i):
    s = (cat['zp'] > ZBINS[i]) & (cat['zp'] <= ZBINS[i + 1])
    return {k: v[s] for k, v in cat.items()}

def calibracion(cat):
    """Responsividad R y sesgo multiplicativo medio m (Mandelbaum+2018)."""
    w = cat['w']
    R = 1.0 - np.sum(w * cat['erms']**2) / np.sum(w)
    m = np.sum(w * cat['m']) / np.sum(w)
    return R, m

def proyectar(ra, dec, ra0=None, dec0=None):
    """Proyeccion gnomonica (plano tangente). x crece hacia el ESTE (RA),
    y hacia el norte. Devuelve grados."""
    ra = np.where(ra - np.median(ra) > 180, ra - 360, ra)
    if ra0 is None: ra0 = 0.5 * (ra.min() + ra.max())
    if dec0 is None: dec0 = 0.5 * (dec.min() + dec.max())
    a, d, a0, d0 = map(np.deg2rad, (ra, dec, ra0, dec0))
    cosc = np.sin(d0)*np.sin(d) + np.cos(d0)*np.cos(d)*np.cos(a - a0)
    x = np.cos(d)*np.sin(a - a0) / cosc
    y = (np.cos(d0)*np.sin(d) - np.sin(d0)*np.cos(d)*np.cos(a - a0)) / cosc
    return np.rad2deg(x), np.rad2deg(y), ra0, dec0

def grilla_shear(cat, pix_arcmin=0.88, R=None, m=None, extent=None):
    """gamma calibrado por pixel:
       gamma = [ sum(w e)/(2 R sum w) - sum(w c)/sum w ] / (1 + m)"""
    if R is None or m is None:
        R, m = calibracion(cat)
    x, y, ra0, dec0 = proyectar(cat['ra'], cat['dec'])
    p = pix_arcmin / 60.0
    if extent is None:
        extent = (x.min(), x.max(), y.min(), y.max())
    nx = int(np.ceil((extent[1] - extent[0]) / p)); ny = int(np.ceil((extent[3] - extent[2]) / p))
    ix = np.clip(((x - extent[0]) / p).astype(int), 0, nx - 1)
    iy = np.clip(((y - extent[2]) / p).astype(int), 0, ny - 1)
    def H(v):
        return np.bincount(iy * nx + ix, weights=v, minlength=nx * ny).reshape(ny, nx)
    w = cat['w']
    W = H(w); N = H(np.ones_like(w))
    with np.errstate(invalid='ignore', divide='ignore'):
        g1 = (H(w * cat['e1']) / (2 * R * W) - H(w * cat['c1']) / W) / (1 + m)
        g2 = (H(w * cat['e2']) / (2 * R * W) - H(w * cat['c2']) / W) / (1 + m)
    g1[W == 0] = 0; g2[W == 0] = 0
    return dict(g1=g1, g2=g2, W=W, N=N, R=R, m=m, extent=extent, pix=pix_arcmin,
                ra0=ra0, dec0=dec0)

def kaiser_squires(g1, g2, pad=2):
    """kappa_E + i kappa_B = F^-1[ D*(k) gamma(k) ],  D = (k1^2 - k2^2 + 2 i k1 k2)/k^2.
    Aqui eje 0 = y, eje 1 = x."""
    ny, nx = g1.shape
    Py, Px = pad * ny, pad * nx
    G = np.zeros((Py, Px), complex); G[:ny, :nx] = g1 + 1j * g2
    k1 = np.fft.fftfreq(Px)[None, :]; k2 = np.fft.fftfreq(Py)[:, None]
    k2tot = k1**2 + k2**2; k2tot[0, 0] = 1.0
    D = (k1**2 - k2**2 + 2j * k1 * k2) / k2tot
    kap = np.fft.ifft2(np.conj(D) * np.fft.fft2(G))[:ny, :nx]
    return kap.real, kap.imag

def suavizar(mapa, peso, sigma_pix):
    """Suavizado gaussiano normalizado (evita bordes/huecos sesgados)."""
    from scipy.ndimage import gaussian_filter
    num = gaussian_filter(mapa * (peso > 0), sigma_pix, mode='constant')
    den = gaussian_filter((peso > 0).astype(float), sigma_pix, mode='constant')
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0.5, num / den, np.nan)
