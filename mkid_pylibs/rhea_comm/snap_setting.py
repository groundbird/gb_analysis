#!/usr/bin/env python3
# coding: utf-8
'''Description of snapshot core handler.'''

from raw_setting import RawSetting

class SnapSetting(RawSetting):
    '''Snapshot handler.'''
    def __init__(self, rbcp_ins, verbose = True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'SNP')
        self.src_list = {'DDS'   : 0,
                         'DAC'   : 1,
                         'ADC'   : 2,
                         'IQ'    : 3,
                         'TEST'  : 4,
                         'PHASE' : 5}
        self.src_list_inv = {str(v):k for k, v in self.src_list.items()}

    def get_snapstatus(self):
        return self._read(0x30000000)

    def snap_on(self, on_off = True):
        data = 1 if on_off else 0
        self._write(0x30000000, data)
        return

    def snap_off(self): return self.snap_on(on_off = False)

    def time_reset(self): self._write(0x30000001, 1)
    def fifo_reset(self): self._write(0x30000002, 1)

    def set_src(self, src, ch = None):
        if not src in self.src_list: return
        self._write(0x31000000, self.src_list[src])
        if type(ch) == int: self._write(0x31000001, ch)
        return

    def get_src(self):
        src_no = str(self._read(0x31000000))
        if not src_no in self.src_list_inv: return "None"
        return self.src_list_inv[src_no]

    def get_ch(self):
        return self._read(0x31000001)

    pass
