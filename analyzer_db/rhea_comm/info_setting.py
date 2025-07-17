#!/usr/bin/env python3
# coding: utf-8
'''Retrieve information stored in the firmware.'''

from raw_setting import RawSetting

INFO_FIRM_VER = 0x00000000
INFO_MAX_CH   = 0x00000010
INFO_EN_SNAP  = 0x00000011
INFO_TRIG_CH  = 0x00000012

class InfoSetting(RawSetting):
    '''Class to retrieve information stored in the firmware.'''
    def __init__(self, rbcp_ins, verbose = True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'INF')

    @property
    def version(self):
        '''Firmware version.

        Returns
        -------
        version : int
            Firmware version.
        '''
        return self._read4(INFO_FIRM_VER)

    @property
    def max_ch(self):
        '''Number of DDS channels.

        Returns
        -------
        max_ch : int
            Number of DDS channels.
        '''
        return self._read(INFO_MAX_CH)

    @property
    def en_snap(self):
        '''Snapshot enabled or not

        Returns
        -------
        en_snap : bool
            True if snapshot enabled.
        '''
        return bool(self._read(INFO_EN_SNAP))

    @property
    def trig_ch(self):
        '''Number of trigger channel.

        Returns
        -------
        trig_ch : int
            Number of trigger enabled channels.
        '''
        return self._read(INFO_TRIG_CH)
