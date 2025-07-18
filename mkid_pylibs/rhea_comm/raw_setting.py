#!/usr/bin/env python3
# coding: utf-8
'''Contains the definition for `RawSetting` class.
'''

class RawSetting:
    '''Generic class for firmware component instances'''
    def __init__(self, rbcp_inst, verbose=True, label='RAW'):
        self._label = label
        self._rbcp_inst = rbcp_inst
        self._verbose = verbose

    def _write_n(self, length, addr, data, byteorder='big', signed=False):
        data_byte = (data).to_bytes(length, byteorder=byteorder, signed=signed)
        self._rbcp_inst.write(addr, data_byte)

        if self._verbose:
            for _i, _d in enumerate(data_byte):
                print(f'{self._label}: write0x{addr+_i:08x} <== 0x{_d:02x} ({_d:3d})')

    def _read_n(self, length, addr, byteorder='big', signed=False):
        data_byte  = self._rbcp_inst.read(addr, length)
        data = int.from_bytes(data_byte, signed=signed, byteorder=byteorder)

        if self._verbose:
            for _i, _d in enumerate(data_byte):
                print(f'{self._label}: write0x{addr+_i:08x} <== 0x{_d:02x} ({_d:3d})')

        return data

    def _write(self, addr, data):
        self._write_n(1, addr, data)

    def _write2(self, addr, data, byteorder='big'):
        self._write_n(2, addr, data, byteorder)

    def _write4(self, addr, data, byteorder='big'):
        self._write_n(4, addr, data, byteorder)

    def _read(self, addr):
        return self._read_n(1, addr)

    def _read2(self, addr, byteorder='big'):
        return self._read_n(2, addr, byteorder)

    def _read4(self, addr, byteorder='big'):
        return self._read_n(4, addr, byteorder)

    @property
    def verbose(self):
        '''Verbosity

        Returns
        -------
        verbose : bool
            True if verbose
        '''
        return self._verbose

    @verbose.setter
    def verbose(self, val):
        assert isinstance(val, bool)
        self._verbose = val

    @property
    def label(self):
        '''Label for the setting.

        Returns
        -------
        label :str
            Label.
        '''
        return self._label
