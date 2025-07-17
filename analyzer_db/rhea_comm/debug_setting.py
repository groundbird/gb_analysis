#!/usr/bin/env python3
# coding: utf-8
'''RBCP debug function
'''

from raw_setting import RawSetting

class DebugSetting(RawSetting):
    '''Class to handle RBCP debugger.'''
    def __init__(self, rbcp_ins, verbose=True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'DBG')

    def probe(self, addr):
        '''Get probe value.'''
        return self._read(0xf1000000 + addr)

    def drive(self, addr, data):
        '''Set drive value.'''
        self._write(0xf2000000 + addr, data)

    def drive_read(self, addr):
        '''Read drive value.'''
        return self._read(0xf2000000 + addr)
