#!/usr/bin/env python3
'''Tone configuration.'''
import numpy as np


class ToneError(Exception):
    '''Error raised in tone configuration.'''


class ToneConf:
    '''Tone configuration.

    Parameters
    ----------
    freq_if_megahz : list of float
        List of tone frequencies in Hz.
    phase : list of float
        Initial tone phases.
    amp : list of float
        Tone amplitudes.
    '''
    def __init__(self, max_ch, freq_if_megahz, phases=None, amps=None, power=1):
        self.freq_if = [f_mega*1e6 for f_mega in freq_if_megahz]
        self._max_ch = max_ch

        if phases is None:
            self.phases = [0.]*self.n_tone
        else:
            self.phases = phases

        if amps is None:
            self.amps = [1.]*self.n_tone
        else:
            self.amps = amps

        if len(self.phases) != self.n_tone:
            raise ToneError('Phase length mismtach.')

        if len(self.amps) != self.n_tone:
            raise ToneError('Amp length mismtach.')

        self.power = power

    @property
    def n_tone(self):
        '''Number of tones.

        Returns
        -------
        n_tone : int
            Number of tones.
        '''
        return len(self.freq_if)

    @property
    def _num_list(self):
        return np.floor(self._max_ch / self.n_tone + 0.1) if self.power < 0 else self.power

    def _mult(self, target):
        if self.power < 0:
            return target * self._num_list
        else:
            tmplist = target * self._num_list
            tmplist += [0.]*(self._max_ch - len(tmplist))
            return tmplist

    @property
    def freq_mult(self):
        '''Frequency list multiplied by `power`'''
        return self._mult(self.freq_if)

    @property
    def amp_mult(self):
        '''Amplitude list multiplied by `power`'''
        return self._mult(self.amps)

    @property
    def phase_mult(self):
        '''Phase list multiplied by `power`'''
        return self._mult(self.phases)

    @property
    def freq_hz_int(self):
        '''Frequency in integer.'''
        return [int(np.floor(freq + 0.5)) for freq in self.freq_if]

    @property
    def freq_repr(self):
        '''Byte representation of frequencies.'''
        brepr = b''
        freq_packs = [int(freq_hz).to_bytes(7, byteorder='big', signed=True) \
                      for freq_hz in self.freq_hz_int]
        for freq_pack in freq_packs:
            brepr += freq_pack * 2

        return brepr

    @property
    def freq_if_megahz(self):
        '''Frequency in MHz.'''
        return [freq/1e6 for freq in self.freq_if]

    def __add__(self, freq_megahz):
        freq = [freq_curr + freq_megahz for freq_curr in self.freq_if_megahz]

        return ToneConf(self._max_ch, freq,
                        phases=self.phases, amps=self.amps, power=self.power)

    def __sub__(self, freq_megahz):
        return self.__add__(-freq_megahz)
