"""
-----------------------------------------------------------------------------------------
ncsf_utils.py
-----------------------------------------------------------------------------------------
Helpers for the nCSF fit with prfmodel (>= 0.2.0):
- stimulus sequence from the events file
- CSF tuning with a fitted left-branch width (width_l), as in the prfpy_csenf fits
- CSF derivatives (log-CSF area, cut-off spatial frequency) without prfpy_csenf
-----------------------------------------------------------------------------------------
"""
import numpy as np
import pandas as pd


def ncsf_stimulus_levels(sf_minFreq, sf_maxFreq, sf_filtNum, minCont, maxCont, contNum):
    """
    Spatial frequency and contrast values of the nCSF task (same as the previous ncsf_fit.py).

    The events file codes the levels with integers starting at 1; the last code is the
    blank (spatial frequency 0, contrast 0).

    Returns
    -------
    sf_levels : np.ndarray, spatial frequencies (cycles/deg), last value 0 (blank)
    con_levels : np.ndarray, Michelson contrasts (%), last value 0 (blank)
    """
    sf_levels = np.concatenate([np.round(np.logspace(np.log10(sf_minFreq), np.log10(sf_maxFreq), sf_filtNum), 2), [0]])
    con_levels = np.concatenate([np.logspace(np.log10(minCont), np.log10(maxCont), contNum), [0]])
    return sf_levels, con_levels


def ncsf_sequence_from_events(events_fn, sf_levels, con_levels):
    """
    Read the per-TR spatial frequency and contrast sequence of an nCSF run.

    Parameters
    ----------
    events_fn : str, path of the *_events.tsv file (one row per TR, columns
                'spatial_frequency' and 'michelson_contrast' with integer codes from 1)
    sf_levels, con_levels : np.ndarray, see ncsf_stimulus_levels

    Returns
    -------
    sf_seq : np.ndarray (n TRs,), spatial frequency in cycles/deg (0 = blank)
    con_seq : np.ndarray (n TRs,), Michelson contrast in % (0 = blank), as prfmodel expects
    """
    events_df = pd.read_table(events_fn, sep="\t")
    sf_map = {i + 1: v for i, v in enumerate(sf_levels)}
    con_map = {i + 1: v for i, v in enumerate(con_levels)}
    sf_seq = events_df['spatial_frequency'].map(sf_map).to_numpy(dtype=float, copy=True)
    con_seq = events_df['michelson_contrast'].map(con_map).to_numpy(dtype=float, copy=True)
    if np.isnan(sf_seq).any() or np.isnan(con_seq).any():
        raise ValueError(f"Unknown spatial frequency / contrast code in {events_fn}")
    # a blank has no grating: zero both, so that prfmodel predicts exactly 0
    blank = (sf_seq == 0) | (con_seq == 0)
    sf_seq[blank], con_seq[blank] = 0.0, 0.0
    return sf_seq, con_seq


def make_csf_tuning_free_width_l():
    """
    CSF tuning (prfmodel CSFTuning) where width_l is a fitted parameter instead of a constant.

    Imported lazily so that this file can be used without keras/prfmodel.
    """
    import math
    from keras import ops
    from prfmodel.models.csf import CSFTuning
    from prfmodel.models.csf._csf import _log10_contrast_sensitivity, _contrast_response

    class CSFTuningFreeWidthL(CSFTuning):
        """CSFTuning with parameters cs_peak, sf_peak, width_l, width_r, slope_crf."""

        @property
        def parameter_names(self):
            return ["cs_peak", "sf_peak", "width_l", "width_r", "slope_crf"]

        def call(self, stimulus, parameters):
            dtype = parameters.dtype
            sf = ops.convert_to_tensor(stimulus.sf, dtype=dtype)
            is_grating = ops.expand_dims(sf > 0, 0)
            log10_sensitivity = _log10_contrast_sensitivity(
                sf=ops.where(sf > 0, sf, ops.ones_like(sf)),
                cs_peak=parameters[["cs_peak"]],
                sf_peak=parameters[["sf_peak"]],
                width_l=parameters[["width_l"]],
                width_r=parameters[["width_r"]],
                dtype=dtype)
            response = _contrast_response(
                contrast=stimulus.contrast,
                log_sensitivity=math.log(10.0) * log10_sensitivity,
                slope_crf=parameters[["slope_crf"]],
                dtype=dtype,
                eps=self.eps)
            return ops.where(is_grating, response, ops.zeros_like(response))

    return CSFTuningFreeWidthL()


def log10_csf(sf, cs_peak, sf_peak, width_l, width_r):
    """
    Asymmetric log-parabolic CSF (log10 sensitivity), same formula as prfmodel/prfpy_csenf.

    sf : (n_sf,) ; parameters : (n_vertices,) or scalars -> returns (n_vertices, n_sf)
    """
    sf = np.asarray(sf, dtype=float)[None, :]
    cs_peak, sf_peak, width_l, width_r = [np.atleast_1d(np.asarray(p, dtype=float))[:, None]
                                          for p in (cs_peak, sf_peak, width_l, width_r)]
    log_diff_sq = (np.log10(sf) - np.log10(sf_peak)) ** 2
    width = np.where(sf < sf_peak, width_l, width_r)
    return np.log10(cs_peak) - log_diff_sq * width ** 2


def ncsf_aulcsf(cs_peak, sf_peak, width_l, width_r, sf_min, sf_max, n_sf=100, normalize=False, cs_max=None):
    """
    Area under the log contrast sensitivity function (AULCSF).

    Integral over log10(SF) in [sf_min, sf_max] of max(log10 CSF, 0) (sensitivity < 1 counts
    as 0), computed on n_sf log-spaced spatial frequencies.

    normalize=True divides by the area of a flat CSF at cs_max over the same SF range
    (log10(cs_max) * (log10 sf_max - log10 sf_min)), so that values lie in [0, 1].
    NB: prfpy_csenf's ncsf_calculate_aulcsf was evaluated on the 6 stimulus SFs only; values are
    therefore not strictly identical to the previous fits.
    """
    sfs = np.logspace(np.log10(sf_min), np.log10(sf_max), n_sf)
    logcsf = np.clip(log10_csf(sfs, cs_peak, sf_peak, width_l, width_r), 0, None)
    auc = np.trapezoid(logcsf, x=np.log10(sfs), axis=1)
    if normalize:
        if cs_max is None:
            raise ValueError("cs_max is required to normalize the AULCSF")
        auc = auc / (np.log10(cs_max) * (np.log10(sf_max) - np.log10(sf_min)))
    return auc


def ncsf_sfmax(cs_peak, sf_peak, width_r, max_sfmax=50):
    """
    Cut-off (acuity) spatial frequency: SF above the peak where the sensitivity falls to 1
    (log10 CSF = 0), i.e. log10(SFmax) = log10(sf_peak) + sqrt(log10(cs_peak)) / width_r.
    Clipped at max_sfmax (same default as prfpy_csenf); vertices with cs_peak <= 1 get sf_peak.
    """
    cs_peak, sf_peak, width_r = [np.asarray(p, dtype=float) for p in (cs_peak, sf_peak, width_r)]
    log_cs = np.clip(np.log10(cs_peak), 0, None)
    sfmax = 10 ** (np.log10(sf_peak) + np.sqrt(log_cs) / width_r)
    return np.clip(sfmax, None, max_sfmax)
