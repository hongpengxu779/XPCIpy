# -*- coding: utf-8 -*-
"""
Phase retrieval algorithms for Propagation-Based Imaging (PBI / Inline).

Supported methods
-----------------
1. **Paganin** (single-distance, homogeneous-material assumption)
   – D. Paganin et al., J. Microsc. 206, 33 (2002)
2. **TIE multi-distance** (Transport of Intensity Equation)
   – generalised Fourier-space solver using ≥2 defocused intensity images
3. **CTF multi-distance** (linearised Contrast Transfer Function)
   – weak-object approximation; better for mixed-phase/absorption objects
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift, fftfreq


# ---------------------------------------------------------------------------
#  Utility helpers
# ---------------------------------------------------------------------------

def _wavelength_um(energy_keV: float) -> float:
    """Return photon wavelength in µm given energy in keV."""
    return 1.23984193 / (energy_keV * 1000.0)


def _freq_grid(ny: int, nx: int, pixel_size_um: float):
    """Return 2-D squared-frequency grid |q|² in µm⁻²."""
    fy = fftfreq(ny, d=pixel_size_um)
    fx = fftfreq(nx, d=pixel_size_um)
    FX, FY = np.meshgrid(fx, fy)
    return FX**2 + FY**2


def _effective_pixel_size(pixel_size_um: float, M: float) -> float:
    """Pixel size in the object plane (µm)."""
    return pixel_size_um / M


def _magnification(dso_cm: float, dod_cm: float) -> float:
    """Geometric magnification for cone beam."""
    if dso_cm <= 0:
        return 1.0
    return (dso_cm + dod_cm) / dso_cm


def _effective_distance_cm(dod_cm: float, M: float) -> float:
    """Effective (reduced) propagation distance for cone beam: z_eff = DOD / M."""
    return dod_cm / M


def normalize_image(raw, flat=None, dark=None, eps=1e-12):
    """
    Flat/dark-field normalisation.

    Parameters
    ----------
    raw : ndarray (float)   – raw intensity image (or stack)
    flat : ndarray or None  – flat-field (open beam)
    dark : ndarray or None  – dark-field (beam off)

    Returns
    -------
    norm : ndarray – normalised intensity  I/I₀
    """
    raw = np.asarray(raw, dtype=np.float64)
    if dark is not None:
        dark = np.asarray(dark, dtype=np.float64)
        raw = raw - dark
    if flat is not None:
        flat = np.asarray(flat, dtype=np.float64)
        if dark is not None:
            flat = flat - dark
        flat = np.where(np.abs(flat) < eps, eps, flat)
        return raw / flat
    return raw


# ---------------------------------------------------------------------------
#  1.  Paganin single-distance retrieval
# ---------------------------------------------------------------------------

def paganin_single_distance(
    image_norm,
    energy_keV: float,
    pixel_size_um: float,
    dod_cm: float,
    delta_beta_ratio: float = 1000.0,
    dso_cm: float = 0.0,
    pad: int = 0,
    regularisation: float = 0.0,
):
    """
    Paganin phase retrieval for a single propagation distance.

    Assumes a *homogeneous* object so that δ(r)/β(r) ≈ const.

    Parameters
    ----------
    image_norm : 2-D ndarray
        Flat/dark corrected **normalised** intensity  I(x,y) / I₀.
    energy_keV : float
        Photon energy in keV.
    pixel_size_um : float
        Detector pixel size in µm.
    dod_cm : float
        Object-to-detector distance in cm.
    delta_beta_ratio : float
        Material δ/β ratio (default 1000; typical for soft tissue at diagnostic energies).
    dso_cm : float
        Source-to-object distance in cm (0 → parallel beam).
    pad : int
        Number of edge-padding pixels (0 = no padding).
    regularisation : float
        Small additive constant in the denominator to avoid division-by-zero (Tikhonov).

    Returns
    -------
    projected_thickness : 2-D ndarray
        μ·T(x,y) = −ln(…); proportional to projected thickness.
    phase : 2-D ndarray
        φ(x,y) = −(2π/λ)·δ·T  (unwrapped, absolute phase shift).
    """
    img = np.asarray(image_norm, dtype=np.float64)
    ny, nx = img.shape

    # Geometry
    M = _magnification(dso_cm, dod_cm) if dso_cm > 0 else 1.0
    pix_obj = _effective_pixel_size(pixel_size_um, M)
    z_eff_cm = _effective_distance_cm(dod_cm, M) if dso_cm > 0 else dod_cm
    z_eff_um = z_eff_cm * 1e4  # cm → µm

    lam = _wavelength_um(energy_keV)

    # Optional padding
    if pad > 0:
        img = np.pad(img, pad, mode='edge')

    Ny, Nx = img.shape
    q2 = _freq_grid(Ny, Nx, pix_obj)  # µm⁻²

    # Paganin filter  (Fourier domain)
    #  -ln{ F⁻¹[ F[I_norm] / (1 + π·λ·z·(δ/β)·|q|²) ] }
    denom = 1.0 + np.pi * lam * z_eff_um * delta_beta_ratio * (4.0 * np.pi**2 * q2)
    denom += regularisation

    F_img = fft2(img)
    filtered = np.real(ifft2(F_img / denom))
    filtered = np.clip(filtered, 1e-30, None)
    projected_thickness = -np.log(filtered)            # µ·T

    # phase = -(δ/β) · µ·T  (sign convention: negative phase for delay)
    phase = -delta_beta_ratio * projected_thickness

    if pad > 0:
        projected_thickness = projected_thickness[pad:-pad, pad:-pad]
        phase = phase[pad:-pad, pad:-pad]

    return projected_thickness, phase


# ---------------------------------------------------------------------------
#  2.  TIE multi-distance retrieval
# ---------------------------------------------------------------------------

def tie_multi_distance(
    images_norm,
    distances_cm,
    energy_keV: float,
    pixel_size_um: float,
    dso_cm: float = 0.0,
    alpha: float = 1e-3,
    pad: int = 0,
):
    """
    Multi-distance TIE (Transport of Intensity Equation) phase retrieval.

    Uses a least-squares fit of dI/dz in Fourier space and inverts the
    Laplacian to recover the phase.

    Parameters
    ----------
    images_norm : list of 2-D ndarrays   (len ≥ 2)
        Normalised intensity images  I_k / I₀  at different DODs.
    distances_cm : list / array of float  (same length as images_norm)
        Corresponding effective propagation distances (DOD) in cm.
    energy_keV : float
    pixel_size_um : float      – detector pixel size in µm
    dso_cm : float             – source-to-object distance (0 → parallel beam)
    alpha : float              – Tikhonov regularisation parameter
    pad : int                  – edge-padding pixels

    Returns
    -------
    phase : 2-D ndarray
        Retrieved phase map  φ(x,y).
    absorption : 2-D ndarray
        Estimated absorption  −ln(I₀) (average of input images).
    dIdz : 2-D ndarray
        Estimated axial intensity derivative (for diagnostics).
    """
    N = len(images_norm)
    if N < 2:
        raise ValueError("TIE multi-distance requires at least 2 images at different distances.")

    imgs = [np.asarray(im, dtype=np.float64) for im in images_norm]
    ny, nx = imgs[0].shape

    # Compute effective distances and pixel
    M_list = []
    zeff_list = []
    for d in distances_cm:
        M = _magnification(dso_cm, d) if dso_cm > 0 else 1.0
        M_list.append(M)
        zeff_list.append(_effective_distance_cm(d, M) if dso_cm > 0 else d)

    # For the TIE we work in the *object plane*; use average M for pixel
    M_avg = float(np.mean(M_list))
    pix_obj = _effective_pixel_size(pixel_size_um, M_avg)

    # Rescale images to object plane intensity (divide by M² for cone beam)
    for k in range(N):
        if M_list[k] != 1.0:
            imgs[k] = imgs[k] / (M_list[k] ** 2)

    zeff_um = np.array(zeff_list) * 1e4  # cm → µm
    lam = _wavelength_um(energy_keV)
    k_wave = 2.0 * np.pi / lam  # wavenumber in µm⁻¹

    # Padding
    if pad > 0:
        imgs = [np.pad(im, pad, mode='edge') for im in imgs]

    Ny, Nx = imgs[0].shape

    # ---------- Estimate dI/dz via least-squares linear fit ----------
    # I_k ≈ I₀ + (dI/dz)·z_k   →  fit slope for each pixel
    z = zeff_um - zeff_um.mean()
    stack = np.stack(imgs, axis=0)                      # (N, Ny, Nx)
    I_mean = stack.mean(axis=0)                         # I₀ estimate

    # slope = Σ (z_k - z̄)(I_k - Ī) / Σ (z_k - z̄)²
    z_centered = z[:, None, None]                       # (N,1,1)
    numerator = np.sum(z_centered * (stack - I_mean[None, ...]), axis=0)
    denominator_z = np.sum(z**2)
    dIdz = numerator / (denominator_z + 1e-30)          # (Ny, Nx)

    # ---------- Solve TIE:  ∇²φ = -(k / I₀) · dI/dz ----------
    q2 = _freq_grid(Ny, Nx, pix_obj)                   # µm⁻²
    laplacian_op = -(2.0 * np.pi) ** 2 * q2             # = -(2π)²|q|²

    rhs = -(k_wave / np.where(np.abs(I_mean) < 1e-30, 1e-30, I_mean)) * dIdz

    F_rhs = fft2(rhs)
    # Regularised inversion of ∇²
    inv_lap = laplacian_op / (laplacian_op**2 + alpha)
    # Set DC to zero (phase offset is arbitrary)
    inv_lap[0, 0] = 0.0
    F_phase = F_rhs * inv_lap
    phase = np.real(ifft2(F_phase))

    absorption = -np.log(np.clip(I_mean, 1e-30, None))

    if pad > 0:
        phase = phase[pad:-pad, pad:-pad]
        absorption = absorption[pad:-pad, pad:-pad]
        dIdz = dIdz[pad:-pad, pad:-pad]

    return phase, absorption, dIdz


# ---------------------------------------------------------------------------
#  3.  CTF multi-distance retrieval  (weak-object approximation)
# ---------------------------------------------------------------------------

def ctf_multi_distance(
    images_norm,
    distances_cm,
    energy_keV: float,
    pixel_size_um: float,
    dso_cm: float = 0.0,
    alpha_phase: float = 1e-2,
    alpha_abs: float = 1e-2,
    pad: int = 0,
):
    """
    Multi-distance CTF (Contrast Transfer Function) phase retrieval.

    Linearised weak-object model (Cloetens et al., 1999):
        I_k(q)/I₀ − 1 ≈ 2·sin(χ_k)·φ(q) − 2·cos(χ_k)·μ(q)
    with  χ_k = π·λ·z_k·|q|².

    Solved via regularised least-squares over multiple distances.

    Parameters
    ----------
    images_norm : list of 2-D ndarrays (len ≥ 2)
    distances_cm : list of float
    energy_keV, pixel_size_um, dso_cm : same as TIE
    alpha_phase, alpha_abs : regularisation for phase and absorption channels
    pad : int

    Returns
    -------
    phase : 2-D ndarray      – retrieved phase φ(x,y)
    absorption : 2-D ndarray  – retrieved absorption μ(x,y)
    """
    N = len(images_norm)
    if N < 2:
        raise ValueError("CTF multi-distance requires at least 2 images.")

    imgs = [np.asarray(im, dtype=np.float64) for im in images_norm]
    ny, nx = imgs[0].shape

    M_list = []
    zeff_list = []
    for d in distances_cm:
        M = _magnification(dso_cm, d) if dso_cm > 0 else 1.0
        M_list.append(M)
        zeff_list.append(_effective_distance_cm(d, M) if dso_cm > 0 else d)

    M_avg = float(np.mean(M_list))
    pix_obj = _effective_pixel_size(pixel_size_um, M_avg)
    zeff_um = np.array(zeff_list) * 1e4
    lam = _wavelength_um(energy_keV)

    for k in range(N):
        if M_list[k] != 1.0:
            imgs[k] = imgs[k] / (M_list[k] ** 2)

    if pad > 0:
        imgs = [np.pad(im, pad, mode='edge') for im in imgs]

    Ny, Nx = imgs[0].shape
    q2 = _freq_grid(Ny, Nx, pix_obj)

    # Build normal-equation accumulators  (2×2 system per frequency pixel)
    A11 = np.zeros((Ny, Nx), dtype=np.float64)
    A12 = np.zeros((Ny, Nx), dtype=np.float64)
    A22 = np.zeros((Ny, Nx), dtype=np.float64)
    b1 = np.zeros((Ny, Nx), dtype=np.complex128)
    b2 = np.zeros((Ny, Nx), dtype=np.complex128)

    for k in range(N):
        chi_k = np.pi * lam * zeff_um[k] * (4.0 * np.pi**2 * q2)
        sin_chi = 2.0 * np.sin(chi_k)
        cos_chi = 2.0 * np.cos(chi_k)

        F_contrast = fft2(imgs[k] - 1.0)  # F[I/I₀ - 1]

        A11 += sin_chi ** 2
        A12 += -sin_chi * cos_chi
        A22 += cos_chi ** 2
        b1 += sin_chi * F_contrast
        b2 += -cos_chi * F_contrast

    # Tikhonov regularisation
    A11 += alpha_phase
    A22 += alpha_abs

    det = A11 * A22 - A12 * A12
    det = np.where(np.abs(det) < 1e-30, 1e-30, det)

    F_phase = (A22 * b1 - A12 * b2) / det
    F_abs = (A11 * b2 - A12 * b1) / det

    phase = np.real(ifft2(F_phase))
    absorption = np.real(ifft2(F_abs))

    if pad > 0:
        phase = phase[pad:-pad, pad:-pad]
        absorption = absorption[pad:-pad, pad:-pad]

    return phase, absorption
