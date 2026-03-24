"""
auto_fitconfig.py
-----------------
Module for automatically determining fit configuration from sweep data.
Used as an alternative to get_fitconfig(chip) inside _perform_fitting.

Return format is identical to get_fitconfig:
  rangeind, freqranges, twokidind, twokidfitter,
  initind, fitinit, depind, dep, skipind, guessskip
"""

import numpy as np
import lmfit
from scipy.signal import find_peaks as scipy_find_peaks
from scipy.ndimage import uniform_filter1d


def _has_double_dip(ikr, prominence_threshold=0.05):
    """
    Detect double dips from sweep data.

    Returns
    -------
    is_double : bool
    dip_freqs : np.ndarray  Detected dip frequencies (sorted)
    """
    amp = ikr.swp.amplitude
    filtered_amp = uniform_filter1d(amp, size=20)
    inv_amp = -filtered_amp
    prominence = prominence_threshold * (np.max(amp) - np.min(amp))
    peaks, _ = scipy_find_peaks(inv_amp, prominence=prominence)

    if len(peaks) == 2:
        peak_amps = amp[peaks]
        sorted_indices = np.argsort(peak_amps)
        peaks = peaks[sorted_indices[:2]]
        dip_freqs = np.sort(ikr.swp.f[peaks])
        return True, dip_freqs
    else:
        dip_freqs = ikr.swp.f[peaks] if len(peaks) > 0 else np.array([ikr.tod.f])
        return False, dip_freqs


def get_auto_fitconfig(kr_list):
    """
    Automatically generate fitconfig from a list of KID analysis objects.
    Returns the same tuple format as get_fitconfig(chip),
    so _perform_fitting can be used without any modification.

    Parameters
    ----------
    kr_list : list of PSDAnalyzer
        List of klib.kidana_psd.PSDAnalyzer objects

    Returns
    -------
    Tuple in the same format as get_fitconfig:
      rangeind, freqranges, twokidind, twokidfitter,
      initind, fitinit, depind, dep, skipind, guessskip
    """
    twokidind    = []
    twokidfitter = []
    initind      = []
    fitinit      = []

    for i, ikr in enumerate(kr_list):
        double, dip_freqs = _has_double_dip(ikr)
        initind.append(i)

        if double:
            center_freq = np.mean(dip_freqs)
            fitter = 'gaolinbg2f' if ikr.tod.f < center_freq else 'gaolinbg2l'
            twokidind.append(i)
            twokidfitter.append(fitter)

            par_val = {
                'fr1':   lmfit.Parameter('fr1',   value=dip_freqs[0],
                                         min=dip_freqs[0] - 2e5, max=dip_freqs[0] + 2e5),
                'fr2':   lmfit.Parameter('fr2',   value=dip_freqs[1],
                                         min=dip_freqs[1] - 2e5, max=dip_freqs[1] + 2e5),
                'Qr1':   lmfit.Parameter('Qr1',   value=1e4,  min=1e3, max=1e5),
                'Qc1':   lmfit.Parameter('Qc1',   value=2e4,  min=1e3, max=1e5),
                'phi01': lmfit.Parameter('phi01', value=0.0,  min=-np.pi, max=np.pi),
                'Qr2':   lmfit.Parameter('Qr2',   value=1e4,  min=1e3, max=1e5),
                'Qc2':   lmfit.Parameter('Qc2',   value=2e4,  min=1e3, max=1e5),
                'phi02': lmfit.Parameter('phi02', value=0.0,  min=-np.pi, max=np.pi),
                'absa':  float(np.max(ikr.swp.amplitude)),
                'arga':  0.0,
                'c':     0.0,
            }
        else:
            fr_val = dip_freqs[0] if len(dip_freqs) > 0 else ikr.tod.f
            par_val = {
                'fr':   lmfit.Parameter('fr',   value=fr_val,
                                        min=fr_val - 2e5, max=fr_val + 2e5),
                'Qr':   lmfit.Parameter('Qr',   value=1e4, min=1e3, max=1e5),
                'Qc':   lmfit.Parameter('Qc',   value=2e4, min=1e3, max=1e5),
                'phi0': lmfit.Parameter('phi0', value=0.0, min=-np.pi, max=np.pi),
            }

        fitinit.append(par_val)

    # rangeind/freqranges, depind/dep, skipind/guessskip are not auto-detected
    # and return empty lists, falling back to defaults in _perform_fitting
    # ([None, None], dep=3, skip=False)
    return [], [], twokidind, twokidfitter, initind, fitinit, [], [], [], []
