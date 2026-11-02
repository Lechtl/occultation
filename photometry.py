import glob
import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from scipy.ndimage import gaussian_filter, median_filter
from scipy.spatial import cKDTree
from skimage.registration import phase_cross_correlation

GAIN = 0.54  # e-/ADU, measured from the flat frames (photon transfer); the header GAIN keyword is a setting
READ_NOISE = 3.15  # ADU per pixel and frame, measured from the bias frames


def list_sequences(directory="20240715_reduced"):
    """Header information of every reduced file, sorted in time."""
    sequences = []
    for file in sorted(glob.glob(os.path.join(directory, "*.science.fits"))):
        header = fits.getheader(file)
        sequences.append({
            "file": file,
            "name": os.path.basename(file)[9:15],
            "n_frames": header["NAXIS3"] if header["NAXIS"] == 3 else 1,
            "exptime": header["EXPTIME"],
            "cycle_time": header.get("KCT", np.nan),
            "mjd_start": header["MJD-OBS"],
        })
    return sequences


def read_frames(sequence):
    """All frames of a sequence as a cube with shape (n_frames, ny, nx)."""
    data = fits.getdata(sequence["file"]).astype(np.float32)
    return data[np.newaxis] if data.ndim == 2 else data


def frame_times(sequence):
    """MJD at mid-exposure of every frame: MJD-OBS (1 s resolution) + i * cycle time + t_exp / 2."""
    index = np.arange(sequence["n_frames"])
    seconds = index * sequence["cycle_time"] + sequence["exptime"] / 2
    return sequence["mjd_start"] + seconds / 86400


def smooth_background(image):
    """Large scale background of an image (NaNs are treated as the image median)."""
    filled = np.where(np.isfinite(image), image, np.nanmedian(image))
    coarse = median_filter(filled[::8, ::8], size=15)
    return np.kron(coarse, np.ones((8, 8)))[:image.shape[0], :image.shape[1]]


def source_mask(shape, positions, radii):
    """True within radii[i] of each (x, y) in positions."""
    mask = np.zeros(shape, dtype=bool)
    for (x, y), radius in zip(positions, radii):
        x0, x1 = max(int(x - radius), 0), min(int(x + radius) + 1, shape[1])
        y0, y1 = max(int(y - radius), 0), min(int(y + radius) + 1, shape[0])
        if x0 >= x1 or y0 >= y1:
            continue
        yy, xx = np.ogrid[y0:y1, x0:x1]
        mask[y0:y1, x0:x1] |= (xx - x) ** 2 + (yy - y) ** 2 <= radius ** 2
    return mask


def background_surface(image, mask, sigma=25.0):
    """
    Smooth background of an image from the pixels outside `mask` (normalised Gaussian smoothing).
    In short exposures a residual pattern of ~0.4 ADU stays on the detector: over the 2800 pixels of
    the aperture 0.1 ADU is 280 ADU, a quarter of the star in 0.3 s, and an annulus cannot remove it.
    """
    good = np.isfinite(image) & ~mask
    weights = gaussian_filter(good.astype(np.float32), sigma)
    surface = gaussian_filter(np.where(good, image, 0).astype(np.float32), sigma)
    return surface / np.where(weights > 1e-3, weights, np.nan)


def file_offset(image, reference):
    """Integer (dx, dy) shift of an image relative to the reference image, from cross-correlation."""
    def prepare(im):
        im = np.where(np.isfinite(im), im, np.nanmedian(im)) - smooth_background(im)
        return gaussian_filter(np.clip(im, np.percentile(im, 1), None), 4)

    (dy, dx), _, _ = phase_cross_correlation(prepare(reference), prepare(image), normalization=None)
    return -dx, -dy


def centroid(image, x, y, half_box=30, smooth=4):
    """Intensity weighted centroid around (x, y), computed on the smoothed image."""
    ix, iy = int(round(x)), int(round(y))
    box = image[iy - half_box:iy + half_box + 1, ix - half_box:ix + half_box + 1]
    box = gaussian_filter(np.where(np.isfinite(box), box, np.nanmedian(box)), smooth)
    box = np.clip(box - np.percentile(box, 50), 0, None)
    box = np.where(box > 0.4 * box.max(), box, 0)
    yy, xx = np.indices(box.shape)
    total = box.sum()
    return (ix - half_box + (xx * box).sum() / total, iy - half_box + (yy * box).sum() / total)


def frame_shifts(frames, positions, start_shift, half_box=30):
    """
    Shift (dx, dy) of each frame relative to the reference image, the median centroid offset
    of the bright stars at `positions`. Offsets are measured in every frame (tracking + image motion).
    """
    shifts = np.zeros((len(frames), 2))
    for i, frame in enumerate(frames):
        offsets = []
        for x, y in positions:
            cx, cy = centroid(frame, x + start_shift[0], y + start_shift[1], half_box)
            offsets.append((cx - x, cy - y))
        shifts[i] = np.median(offsets, axis=0)
    return shifts


def measure_fwhm(image, positions, half_box=40):
    """Median FWHM [px] of the stars at `positions`, from their second moment."""
    fwhm = []
    for x, y in positions:
        ix, iy = int(round(x)), int(round(y))
        box = image[iy - half_box:iy + half_box + 1, ix - half_box:ix + half_box + 1]
        box = np.where(np.isfinite(box), box, np.nanmedian(box))
        box = gaussian_filter(box, 1.5) - np.median(box)
        box = np.where(box > 0.2 * box.max(), box, 0)
        yy, xx = np.indices(box.shape)
        total = box.sum()
        mx, my = (xx * box).sum() / total, (yy * box).sum() / total
        var = (((xx - mx) ** 2 + (yy - my) ** 2) * box).sum() / total / 2
        fwhm.append(2.355 * np.sqrt(var))
    return float(np.median(fwhm))


def aperture_photometry(image, x, y, radius, annulus=None):
    """
    Flux in a circular aperture. Bad pixels (NaN) inside the aperture are replaced from their neighbours.
    With an annulus (r_in, r_out) its sigma-clipped mean is subtracted, without one the image must
    already be background-subtracted.
    Returns the flux [ADU], the background per pixel [ADU] and the pixels used for the background.
    """
    half = int(annulus[1] if annulus else np.max(radius)) + 8
    ix, iy = int(round(x)), int(round(y))
    if ix - half < 0 or iy - half < 0 or ix + half >= image.shape[1] or iy + half >= image.shape[0]:
        return np.nan, np.nan, 0
    cutout = image[iy - half:iy + half + 1, ix - half:ix + half + 1]
    px, py = x - ix + half, y - iy + half

    background, n_sky = 0.0, 0
    if annulus:
        yy, xx = np.indices(cutout.shape)
        r = np.hypot(xx - px, yy - py)
        sky = sigma_clip(cutout[(r >= annulus[0]) & (r <= annulus[1])], sigma=3, maxiters=3, masked=True, cenfunc='median')
        background = float(sky.mean())
        n_sky = int(np.sum(np.isfinite(sky.compressed())))

    good = np.isfinite(cutout)
    filled = np.where(good, cutout, 0.0)
    weights = gaussian_filter(good.astype(np.float32), 2.0)
    smoothed = gaussian_filter(filled, 2.0) / np.where(weights > 0, weights, 1)
    filled = np.where(good, cutout, smoothed) - background

    radii = np.atleast_1d(radius)
    flux = [CircularAperture((px, py), r=rad).do_photometry(filled, method='exact')[0][0] for rad in radii]
    return (float(flux[0]) if np.ndim(radius) == 0 else np.array(flux)), background, n_sky


def flux_error(flux, background, radius, n_sky, n_frames=1):
    """1 sigma error [ADU] of an aperture flux on a sum of n_frames frames (CCD equation)."""
    area = np.pi * radius ** 2
    pixel_variance = np.clip(background, 0, None) / GAIN + n_frames * READ_NOISE ** 2
    return np.sqrt(np.clip(flux, 0, None) / GAIN + area * pixel_variance + area ** 2 * pixel_variance / max(n_sky, 1))


def stack_frames(frames, shifts):
    """
    Align frames on the stars with integer pixel shifts and add them up, without interpolating.
    A pixel that is NaN in some frames is scaled up from the frames where it is valid.
    Returns the summed image [ADU]; the stars end up at their position in the reference image.
    """
    total = np.zeros_like(frames[0], dtype=np.float32)
    count = np.zeros_like(frames[0], dtype=np.float32)
    for frame, (dx, dy) in zip(frames, shifts):
        aligned = np.roll(frame, (-int(round(dy)), -int(round(dx))), axis=(0, 1))
        good = np.isfinite(aligned)
        total[good] += aligned[good]
        count[good] += 1
    with np.errstate(invalid='ignore', divide='ignore'):
        return total / count * len(frames)


def solve_field(image, catalog, header, n_stars=60):
    """
    Plate solution of an image against a catalog with columns ra, dec. The header WCS is ~2 arcmin
    off, but its scale and orientation only need a small correction: search rotation, scale and shift
    with a 2D histogram of position differences, then fit an affine transformation to the matches.
    Returns the catalog positions in the image (x, y) and the rms of the fit.
    """
    filled = np.where(np.isfinite(image), image, np.nanmedian(image))
    smooth = gaussian_filter(filled - smooth_background(image), 4)
    threshold = 8 * 1.4826 * np.median(np.abs(smooth - np.median(smooth)))
    sources = DAOStarFinder(fwhm=14, threshold=threshold, exclude_border=True)(smooth)
    sources = sources[np.argsort(-sources['flux'])][:n_stars]
    x, y = np.array(sources['xcentroid']), np.array(sources['ycentroid'])
    all_x, all_y = WCS(header, naxis=2).all_world2pix(catalog['ra'], catalog['dec'], 0)
    bright = np.array(catalog['phot_g_mean_mag']) < 18
    gx, gy = all_x[bright], all_y[bright]

    best = (0,)
    for rotation in np.arange(-6, 6.01, 0.5):
        for scale in [0.98, 0.99, 1.0, 1.01, 1.02]:
            c, s = np.cos(np.radians(rotation)), np.sin(np.radians(rotation))
            px = 1024 + scale * (c * (gx - 1024) - s * (gy - 1024))
            py = 1024 + scale * (s * (gx - 1024) + c * (gy - 1024))
            hist, xedges, yedges = np.histogram2d((x[:, None] - px).ravel(), (y[:, None] - py).ravel(),
                                                  bins=[np.arange(-2500, 2500, 8)] * 2)
            i, j = np.unravel_index(hist.argmax(), hist.shape)
            if hist[i, j] > best[0]:
                best = (hist[i, j], px + xedges[i] + 4, py + yedges[j] + 4)

    distance, index = cKDTree(np.c_[best[1], best[2]]).query(np.c_[x, y], distance_upper_bound=15)
    matched = np.isfinite(distance)
    design = np.c_[gx[index[matched]], gy[index[matched]], np.ones(matched.sum())]
    coeff_x = np.linalg.lstsq(design, x[matched], rcond=None)[0]
    coeff_y = np.linalg.lstsq(design, y[matched], rcond=None)[0]
    rms = np.hypot(np.std(design @ coeff_x - x[matched]), np.std(design @ coeff_y - y[matched]))
    all_stars = np.c_[all_x, all_y, np.ones(len(all_x))]
    return all_stars @ coeff_x, all_stars @ coeff_y, rms, int(matched.sum())
