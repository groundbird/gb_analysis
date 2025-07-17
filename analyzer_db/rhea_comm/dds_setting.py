#!/usr/bin/env python3
# coding: utf-8
'''Description of the class to handle DDS core inside FPGA.
'''
from math import pi

from raw_setting import RawSetting
from rhea_pkg import FREQ_CLK_HZ, PHASE_WIDTH, DDS_AMP_BW
from tone_conf import ToneConf

TRIGGER_ENABLE  = 0x40000000
PHASE_RESET     = 0x40000001
DDS_OFFSET      = 0x41000000
DDS_PERIOD      = 0x41010000

DDS_PINC = lambda channel: DDS_OFFSET + (channel << 8) + (0 << 4)
DDS_POFF = lambda channel: DDS_OFFSET + (channel << 8) + (1 << 4)
DDS_AMPI = lambda channel: DDS_OFFSET + (channel << 8) + (2 << 4)


def pinc2freq(pinc):
    '''Convert pinc value to frequency in Hz.

    Parameter
    ---------
    pinc : int
        Phase increment per clock cycle.

    Returns
    -------
    freq : float
        Frequency in Hz.
    '''
    return pinc * float(FREQ_CLK_HZ) / (1 << PHASE_WIDTH)

def freq2pinc(freq):
    '''Convert frequency in Hz to pinc value.

    Parameter
    ---------
    freq : float
        Frequency in Hz.

    Returns
    -------
    pinc : int
        Phase increment per clock cycle.
    '''
    return int(float(freq) / FREQ_CLK_HZ * (1 << PHASE_WIDTH))

def poff2rad(poff):
    '''Convert poff value to radian.

    Parameter
    ---------
    poff : int
        Phase offset in integer.

    Returns
    -------
    phase : float
        Phase offset in radian.
    '''
    return poff * 2 * pi / (1 << PHASE_WIDTH)

def rad2poff(phase):
    '''Convert radian to poff

    Parameter
    ---------
    phase : float
        Phase offset in radian.

    Returns
    -------
    poff : int
        Phase offset in integer.
    '''
    return int(phase / 2 / pi * (1<<32))

def amp2ampi(amp):
    '''Convert float amplitude to amp_int

    Parameter
    ---------
    amp : float
        Amplitude whose maximum is 1.

    Returns
    -------
    amp_int : int
        Integer representation of the given amp.
    '''
    return int(float(amp) * ((1 << DDS_AMP_BW) - 1))

def ampi2amp(amp_int):
    '''Convert amp_int to float amplitude.

    Parameter
    ---------
    amp_int : int
        Integer representation of the given amp.

    Returns
    -------
    amp : float
        Amplitude whose maximum is 1.
    '''
    return float(amp_int) / ((1 << DDS_AMP_BW) - 1)


class DdsSetting(RawSetting):
    '''Class for DDS setting.'''
    def __init__(self, rbcp_ins, max_ch=2, verbose=True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'DDS')
        self.max_ch = max_ch

    def trig_enable(self):
        '''Enable trigger.'''
        self._write(TRIGGER_ENABLE, 0x01)

    def set_periodic_sync(self, do_sync):
        '''Turn on/off the periodic sync.
        The firmware has a functionality to reset DDS phase once per 1 ms.
        This function turns on/off the phase reset.

        Parameter
        ---------
        do_sync : bool
            Do phase reset.
        '''
        self._write(PHASE_RESET, int(do_sync))

    def get_periodic_sync(self):
        '''Get periodic sync configuration

        Returns
        -------
        do_sync : bool
            Do phase reset if True.
        '''
        return bool(self._read(PHASE_RESET))

    def set_pinc(self, channel, pinc, dds_reset=True):
        '''Set PINC (phase incremental) value to the DDS.
        Frequency of the DDS is specified by the phase incremental value per each clock cycle.

        Parameters
        ----------
        channel : int
            DDS channel.
        pinc : int, float
            Phase incremental value.
        dds_reset : bool, optional
            Do DDS reset after configuration. (Default: True)
        '''
        assert 0 <= channel < self.max_ch

        data  = int(pinc)
        data %= (1<<32)
        self._write4(DDS_PINC(channel), data)

        if dds_reset:
            self.trig_enable()

    def get_pinc(self, channel):
        '''Get PINC value.

        Parameter
        ---------
        channel : int
            DDS channel

        Returns
        -------
        pinc : int
            PINC value
        '''
        assert 0 <= channel < self.max_ch

        data = self._read4(DDS_PINC(channel))
        return data

    def set_poff(self, channel, poff, dds_reset=True):
        '''Set phase offset value.

        Parameters
        ----------
        channel : int
            DDS channel
        poff : int, float
            Phase offset value.
        dds_reset : bool
            Do DDS reset after the configuration if True.
        '''
        assert 0 <= channel < self.max_ch

        data  = int(poff)
        data %= (1<<32)
        self._write4(DDS_POFF(channel), data)

        if dds_reset:
            self.trig_enable()


    def get_poff(self, channel):
        '''Get phase offset value.

        Parameter
        ---------
        channel : int
            DDS channel.

        Returns
        -------
        poff : int
            Phase offset value.
        '''
        assert 0 <= channel < self.max_ch
        poff = self._read4(DDS_POFF(channel))

        return poff

    def set_ampi(self, channel, amp_int, dds_reset=True):
        '''Set amplitude of the channel with raw integer value.

        Parameters
        ----------
        channel : int
            DDS channel
        amp_int : int
            Phase offset value.
            The value should be from 0 (min) to 2**DDS_AMP_BW - 1 (max).
        dds_reset : bool
            Do DDS reset after the configuration if True.
        '''

        assert 0 <= channel < self.max_ch
        assert 0 <= amp_int < 2**DDS_AMP_BW

        self._write4(DDS_AMPI(channel), amp_int)

        if dds_reset:
            self.trig_enable()

    def get_ampi(self, channel):
        '''Get amp value in integer.

        Parameter
        ---------
        channel : int
            DDS channel.

        Returns
        -------
        amp_int : int
            Amplitude in integer.
        '''
        assert 0 <= channel < self.max_ch
        amp_int = self._read4(DDS_AMPI(channel))

        return amp_int

    def set_sync_span(self, span):
        '''Set period of DDS phase reset.

        Parameters
        ----------
        span : int
            Period of DDS phase reset.
        '''
        assert isinstance(span, int)
        assert 0 <= span < (1<<32)

        if span == 0:
            self.set_periodic_sync(False)

        self._write4(DDS_PERIOD, span)
        self.set_periodic_sync(True)

    def get_sync_span(self):
        '''Get period of DDS phse reset.

        Returns
        -------
        span : int
            Period of DDS phase reset.
        '''
        return self._read4(DDS_PERIOD)

    def set_freq(self, channel, freq, dds_reset=True):
        '''Configure DDS frequency using frequency in Hz.
        Parameters
        ----------
        channel : int
            DDS channel.
        freq : float
            Frequency in Hz.
        dds_reset : bool, optional
            Do DDS reset after configuration. (Default: True)
        '''

        self.set_pinc(channel, freq2pinc(freq), dds_reset)

    def set_freqs(self, freq_list):
        '''Configure DDS frequencies according to the given list.

        Parameters
        ----------
        freq_list : list of float
            List of frequencies in Hz.
        '''

        assert len(freq_list) == self.max_ch

        for channel, freq in enumerate(freq_list):
            self.set_freq(channel, freq, dds_reset=False)

        self.trig_enable()

    def get_freq(self, channel):
        '''Get DDS frequency.

        Parameter
        ---------
        channel : int
            DDS channel

        Returns
        -------
        freq : float
            Frequency in Hz.
        '''
        pinc = self.get_pinc(channel)
        return pinc2freq(pinc)

    def set_phase(self, channel, phase, dds_reset=True):
        '''Set initial phase to the DDS.

        Parameters
        ----------
        channel : int
            DDS channel
        phase : float
            Phase offset value in radian.
        dds_reset : bool
            Do DDS reset after the configuration if True.
        '''
        self.set_poff(channel, rad2poff(phase), dds_reset)

    def get_phase(self, channel):
        '''Get initial phase of the DDS.

        Parameter
        ---------
        channel : int
            DDS channel

        Returns
        -------
        phase : float
            Phase offset value in radian.
        '''
        poff = self.get_poff(channel)

        return poff2rad(poff)

    def set_amp(self, channel, amp, dds_reset=True):
        '''Set amplitude of the channel with a float value
        whose maximum is normalized to 1.

        Parameters
        ----------
        channel : int
            DDS channel
        amp : float
            Phase offset value.
        dds_reset : bool
            Do DDS reset after the configuration if True.
        '''
        self.set_ampi(channel, amp2ampi(amp), dds_reset)

    def get_amp(self, channel):
        '''Get amp value in float.

        Parameter
        ---------
        channel : int
            DDS channel.

        Returns
        -------
        amp : float
            Float amplitude.
        '''
        amp_int = self.get_ampi(channel)

        return ampi2amp(amp_int)

    def set_amps(self, amp_list):
        '''Configure DDS amplitudes according to the given list.

        Parameters
        ----------
        amp_list : list of float
            List of amplitudes.
        '''
        assert len(amp_list) == self.max_ch

        for channel, amp in enumerate(amp_list):
            self.set_amp(channel, amp, dds_reset=False)

        self.trig_enable()

    def set_phases(self, phase_list):
        '''Configure DDS initial phases according to the given list.

        Parameters
        ----------
        phase_list : list of float
            List of phases in radian.
        '''
        assert len(phase_list) == self.max_ch

        for channel, phase in enumerate(phase_list):
            self.set_phase(channel, phase, dds_reset=False)

        self.trig_enable()

    def reset(self):
        '''Reset all DDSs'''
        for channel in range(self.max_ch):
            self.set_freq(channel, 0., False)
            self.set_phase(channel, 0., False)

        self.set_sync_span(200000) # 1 kHz span
        self.trig_enable()

    def configure(self, tone_conf:ToneConf):
        '''Configure DDS using `ToneConf` class.

        Parameter
        ---------
        tone_conf : ToneConf
            Tone configuration.
        '''
        self.set_freqs(tone_conf.freq_mult)
        self.set_amps(tone_conf.amp_mult)
        self.set_phases(tone_conf.phase_mult)
