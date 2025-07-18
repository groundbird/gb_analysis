#!/usr/bin/env python3
# coding: utf-8
'''Trigger core setting'''

from raw_setting import RawSetting

TRG_STATUS   = 0x70000000
TRG_POSITION = 0x70000010
TRG_THRCOUNT = 0x70000020
TRG_PRE_MAX = 1000
TRG_ENABLE = lambda ch: 0x71000000 + (ch << 8)
TRG_THR_MAX = (1<<55)
TRG_THR_MIN = -(1<<55) + 1
TRG_THR_BASE = 0x71000000


class TrgSetting(RawSetting):
    '''Trigger core setting'''
    def __init__(self, rbcp_ins, max_ch=2, verbose=True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'TRG')
        self.max_ch = max_ch

    def init(self):
        '''Initialize the core with default values.'''
        for i in range(self.max_ch):
            self.set_threshold(i, [0, 0], [0, 0],
                               trg_reset=False, set_en=False)
            self.set_disable(i, trg_reset=False)

        self.set_trig_pos(100, trg_reset = False)
        self.set_thre_count(1, trg_reset = True)

    def reset(self):
        '''Reset the measurement.'''
        self._write(TRG_STATUS, 0)

    def start(self):
        '''Start self-trigger measurement.'''
        self._write(TRG_STATUS, 1)

    def state(self):
        '''Self-trigger measuemrent status

        Returns
        -------
        running : bool
            True if trigger measurement is going on.
            Goes False when the measurement finishes.
        '''
        return bool(self._read(TRG_STATUS))

    def set_trig_pos(self, head_length=100, trg_reset=True):
        '''Set trigger position.

        Parameters
        ----------
        head_length : int
            Length of data points ahead of triggering.
        trg_reset : bool, optional
            Do reset after the configuration (default: True).
        '''

        assert isinstance(head_length, int), 'set_trig_pos: head_length should be int.'
        assert 0 < head_length < TRG_PRE_MAX, 'set_trig_pos: head_length is 0 -- 1000'

        self._write2(TRG_POSITION, head_length)
        if trg_reset:
            self.reset()

    def get_trig_pos(self):
        '''Get trigger position.

        Returns
        -------
        trig_pos : int
            Length of data points ahead of triggering.
        '''
        return self._read2(TRG_POSITION)

    def set_enable(self, channel, enable=True, trg_reset=True):
        '''Select what channels to be used as trigger sources.
        If you enable more than one channels as trigger sources,
        the self-trigger will be issued when more than one channel
        surpasses the threshold value (i.e. OR logic adopted).

        Parameters
        ----------
        channel : int
            Channel number to be used as a trigger source.
        enable : bool, optional
            Enable (disable) if True (False) (default: True).
        trg_reset : bool
            Do reset after the configuration (default: True).
        '''
        assert channel in range(self.max_ch), f'Ch{channel:d} does not exist.'

        self._write(TRG_ENABLE(channel), int(enable))

        if trg_reset:
            self.reset()

    def set_disable(self, channel, trg_reset=True):
        '''Remove specified channel from trigger source list.'''
        self.set_enable(channel, enable=False, trg_reset=trg_reset)

    def get_enable(self, channel):
        '''Tells whether the channel is in trigger source list.

        Parameter
        ---------
        channel : int
            Channel number.
        '''
        assert channel in range(self.max_ch), f'Ch{channel:d} does not exist.'

        return bool(self._read(TRG_ENABLE(channel)))

    def set_threshold(self, channel, i_val_range, q_val_range,
                      trg_reset=True, set_en=True):
        '''Set threshold for self triggering.

        Parameters
        ----------
        channel : int
            Target channel number.
        i_val_range : (int, int)
            Upper and lower threshold for the trigger for I channel.
        q_val_range : (int, int)
            Upper and lower threshold for the trigger for Q channel.
        trg_reset : bool
            Do reset after the configuration (default: True).
        set_en : bool
            Set the channel enable after the configuration (default: True).
        '''
        assert channel in range(self.max_ch), f'Ch{channel:d} does not exist.'

        val_range = [i_val_range[0],
                     q_val_range[0],
                     i_val_range[1],
                     q_val_range[1]]

        for i, thr_val in enumerate(val_range):
            assert isinstance(thr_val, int), 'set_threshold: min/max_val should be int.'
            assert TRG_THR_MIN <= thr_val <= TRG_THR_MAX,\
            'set_threshold: min/max_val should be 56-bit-int.'

            addr  = TRG_THR_BASE
            addr += (channel << 8)
            addr += ((i+1) << 4)

            self._write_n(8, addr, thr_val, signed=True)

        if set_en:
            self.set_enable(channel, trg_reset=False)

        if trg_reset:
            self.reset()

    def get_threshold(self, channel):
        '''Get threshold value.

        Parameters
        ----------
        channel : int
           Target channel number.

        Returns
        -------
        threshold : [[int, int], [int, int]]
            [[i_min, q_min], [i_max, q_max]]
        '''
        assert channel in range(self.max_ch), f'Ch{channel:d} does not exist.'

        ret = [0, 0, 0, 0]
        for i in range(4):
            addr  = TRG_THR_BASE
            addr += (channel << 8)
            addr += ((i+1) << 4)
            
            ret[i] = self._read_n(8, addr, signed=True)

        return [ret[0:2], ret[2:4]]

    def set_thre_count(self, count, trg_reset=True):
        '''Sets the number of points exceeding threashold reqired to issue trigger.
        The trigger is issued if the data exceeds the threshold for the specified number
        of consective points.

        Parameters
        ----------
        count : int
            Number of threshold-exceeding points required for issuing trigger.
        trg_reset : bool
            Do reset after the configuration (default: True).
        '''
        assert isinstance(count, int), 'set_thre_count: count should be int.'
        assert 0 < count < 1000, 'set_thre_count: count should be 0 -- 1000'

        self._write2(0x70000020, count & 0xffff)

        if trg_reset:
            self.reset()

    def get_thre_count(self):
        '''Gets the number of points exceeding threashold reqired to issue trigger.

        Returns
        -------
        count : int
            Number of threshold-exceeding points required for issuing trigger.
        '''
        return self._read2(0x70000020)
