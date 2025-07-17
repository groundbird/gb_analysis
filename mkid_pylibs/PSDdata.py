################################################################
# PSD data
################################################################
from .BaseMultiData import BaseMultiData
from .FreqTimeArray import FreqArray
from .GeneralArray import named_array

AmpPhsArray = named_array("AmpPhsArray", ['amplitude', 'phase'])

class PSDdata(BaseMultiData,FreqArray,AmpPhsArray):
    """
    A class for PSD data: FreqArray + AmpPhsArray.
    freq[Hz] vs Amplitude/Phase

    - **PSDdata(frq, amp, phs)**:
    :param frq: 1-D array of float of frequency [Hz]
    :param amp: 1-D array of amplitude
    :param phs: 1-D array of phase
    """
    def __init__(self, *args, **kws):
        if args:
            if 'info' in kws:
                info = kws['info']
            else:
                info = dict()

            xdata = FreqArray(args[0])
            ydata = AmpPhsArray((args[1], args[2]))
            super().__init__((xdata, ydata), info)

