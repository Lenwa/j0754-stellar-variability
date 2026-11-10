"""Forced density measurements; scan values are not incident photon flux."""
import numpy as np


def crop_origin(x, y, nx, ny, size, clamp=False):
    if size <= 0 or size > min(nx,ny):
        raise ValueError("Invalid crop size")
    x0, y0 = int(round(x))-size//2, int(round(y))-size//2
    if clamp:
        x0, y0 = max(0,min(x0,nx-size)), max(0,min(y0,ny-size))
    if not (0 <= x0 and 0 <= y0 and x0+size <= nx and y0+size <= ny):
        raise ValueError("Requested field is not fully inside scan")
    return x0, y0


def blank_positions(x, y, width, height, count=400, seed=20261108):
    finite = np.isfinite(x+y)
    x, y = x[finite], y[finite]
    rng = np.random.default_rng(seed)
    result = []
    for bx, by in rng.uniform([25,25], [width-25,height-25], size=(20000,2)):
        if not len(x) or np.min((x-bx)**2+(y-by)**2) > 25**2:
            result.append((float(bx), float(by)))
        if len(result) == count:
            break
    return result


def decode_band(raw, header, x0, y0, width, height):
    """Decode complete FITS rows, then crop columns, preserving the sky mapping."""
    if header["BITPIX"] != 16 or header["NAXIS"] != 2:
        raise ValueError("Only two-dimensional 16-bit plate scans are supported")
    nx = header["NAXIS1"]
    data = np.frombuffer(raw, dtype=">i2").reshape(height, nx)
    crop = data[:, x0:x0+width].astype(float)
    crop = crop * header.get("BSCALE", 1) + header.get("BZERO", 0)
    h = header.copy()
    for key in ["BZERO", "BSCALE", "CHECKSUM", "DATASUM"]:
        h.remove(key, ignore_missing=True)
    h["NAXIS1"], h["NAXIS2"] = width, height
    h["CRPIX1"] -= x0
    h["CRPIX2"] -= y0
    h.add_history(f"HTTP row-band crop: x0={x0}, y0={y0}, width={width}, height={height}")
    return crop, h


def aperture(image, x, y, radius=6, polarity=1):
    """Fixed-position aperture with robust local annular background; no clipping."""
    half = 21
    ix, iy = int(round(x)), int(round(y))
    if ix-half < 0 or iy-half < 0 or ix+half >= image.shape[1] or iy+half >= image.shape[0]:
        return None
    patch = image[iy-half:iy+half+1, ix-half:ix+half+1]
    yy, xx = np.indices(patch.shape)
    rr = np.hypot(xx+ix-half-x, yy+iy-half-y)
    sky = patch[(rr >= 12) & (rr <= 20)]
    background = np.median(sky)
    noise = 1.4826*np.median(np.abs(sky-background))
    weights = rr <= radius
    signal = polarity * np.sum(patch[weights]-background)
    return dict(signal=float(signal), background=float(background),
                pixel_noise=float(noise), area=int(weights.sum()))


def fit_calibration(signals, magnitudes):
    """Local empirical quadratic magnitude(log density), with bounded support."""
    signals, magnitudes = np.asarray(signals), np.asarray(magnitudes)
    good = np.isfinite(signals) & np.isfinite(magnitudes) & (signals > 0)
    if good.sum() < 12 or np.ptp(magnitudes[good]) < 1:
        return None
    x, y = np.log10(signals[good]), magnitudes[good]
    keep = np.ones(len(x), dtype=bool)
    for _ in range(4):
        if keep.sum() < 12 or np.ptp(y[keep]) < 1:
            return None
        coefficients = np.polyfit(x[keep], y[keep], 2)
        residual = y-np.polyval(coefficients, x)
        scatter = max(.05, 1.4826*np.median(np.abs(residual[keep]-np.median(residual[keep]))))
        keep = np.abs(residual-np.median(residual[keep])) <= 3*scatter
    if keep.sum() < 12 or np.ptp(y[keep]) < 1:
        return None
    coefficients = np.polyfit(x[keep], y[keep], 2)
    low, high = float(x[keep].min()), float(x[keep].max())
    if np.max(np.polyval(np.polyder(coefficients), [low, high])) >= 0:
        return None
    return dict(coefficients=coefficients.tolist(), log_signal_min=low,
                log_signal_max=high, stars=int(keep.sum()))


def predict_magnitude(model, signal):
    if model is None or not np.isfinite(signal) or signal <= 0:
        return None
    x = np.log10(signal)
    if not model["log_signal_min"] <= x <= model["log_signal_max"]:
        return None
    if model.get("kind") == "isotonic":
        return float(np.interp(x, model["x"], model["y"]))
    return float(np.polyval(model["coefficients"], x))


def fit_isotonic(signals, magnitudes):
    """Secondary monotonic empirical fit; retained validation decides usefulness."""
    signals, magnitudes = np.asarray(signals), np.asarray(magnitudes)
    good = np.isfinite(signals) & np.isfinite(magnitudes) & (signals > 0)
    if good.sum() < 12 or np.ptp(magnitudes[good]) < 1:
        return None
    x, y = np.log10(signals[good]), magnitudes[good]
    ux, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    means = np.array([np.mean(y[inverse == j]) for j in range(len(ux))])
    if len(ux) < 4:
        return None
    # Pool-adjacent-violators for a non-increasing magnitude response.
    blocks = []
    for i, value in enumerate(means):
        blocks.append([i, i+1, float(value), int(counts[i])])
        while len(blocks) > 1 and blocks[-2][2] < blocks[-1][2]:
            a, b = blocks[-2:]
            weight = a[3]+b[3]
            blocks[-2:] = [[a[0], b[1], (a[2]*a[3]+b[2]*b[3])/weight, weight]]
    fitted = np.empty(len(ux))
    for lo, hi, value, _ in blocks:
        fitted[lo:hi] = value
    return dict(kind="isotonic", x=ux.tolist(), y=fitted.tolist(),
                log_signal_min=float(ux[0]), log_signal_max=float(ux[-1]), stars=int(good.sum()))
