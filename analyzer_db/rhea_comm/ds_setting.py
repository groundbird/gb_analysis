#!/usr/bin/env python3
# coding: utf-8
'''Definition for downsampler handler'''

from raw_setting import RawSetting

DS_OFFSET = 0x61000000

class DsSetting(RawSetting):
    '''Downsampler handler'''
    def __init__(self, rbcp_ins, verbose=True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'DS')

    def set_accum(self, accum_num):
        '''Set accumulation number.

        Parameter
        ---------
        accum_num : int
            Accumulation number.
        '''
        assert 10 <= accum_num <= 200000
        self._write4(DS_OFFSET, accum_num)

    def get_accum(self):
        '''Get accumulation number.

        Returns
        -------
        accum_num : int
            Accumulation number.
        '''
        return self._read4(DS_OFFSET)
