#!/usr/bin/env python3
# coding: utf-8
'''IQ reader module.'''
from raw_setting import RawSetting

IQ_STATUS     = 0x50000000
IQ_RESET_TS   = 0x50000001
IQ_FIFO_ERR   = 0x50000002
IQ_READ_WIDTH = 0x50000010


class IQSetting(RawSetting):
    '''Class for IQ reader module in the firmware.'''
    def __init__(self, rbcp_ins, verbose =True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'IQ ')

    @property
    def iq_status(self):
        '''Check whether the IQ stream output is on or not.

        Returns
        -------
        iq_status : bool
            True if the stream is valid.
        '''
        return bool(self._read(IQ_STATUS))

    def iq_on(self, on_off=True):
        '''Turn on/off IQ stream.

        Parameter
        ---------
        on_off : bool, optional
            True to turn on the stream (default: True).
        '''
        data = 1 if on_off else 0
        self._write(IQ_STATUS, data)

    def iq_off(self):
        '''Turn off IQ stream.'''
        return self.iq_on(on_off=False)

    def time_reset(self):
        '''Reset timestamp.'''
        self._write(IQ_RESET_TS, 1)

    @property
    def fifo_error(self):
        '''FIFO error status.

        Returns
        -------
        fifo_error : bool
            True when fifo error was detected.
        '''
        return bool(self._read(IQ_FIFO_ERR))

    def clear_fifo_error(self):
        '''Clear FIFO error status.'''
        self._write(IQ_FIFO_ERR, 0)

    def set_read_width(self, ch_width):
        '''Set how many channels to be read.
        This functionality provides a way to reduce data trafic
        when the number of channels is smaller than maximum.

        Parameter
        ---------
        ch_width : int
            Number of channels to be read.
        '''
        self._write(IQ_READ_WIDTH, ch_width)

    @property
    def read_width(self):
        '''Number of channels to be read.

        Returns
        -------
        ch_width : int
            Number of channels to be read.
        '''
        return self._read(IQ_READ_WIDTH)
