#!/usr/bin/env python3
# coding: utf-8
'''Description on the ADC setting.'''
from enum import IntEnum
from raw_setting import RawSetting

ADC_SWAP_ADDR = 0x1200_0000
ADS_OFFSET = 0x1000_0000

ADS_TEST_PATTERN_MASK = 0b0000111

class ADSWriteRead:
    '''ADS write/read setting.'''
    write = True
    read = False

class ADSChannel(IntEnum):
    '''ADS channel'''
    A = 0
    B = 1

class ADSTestPattern(IntEnum):
    '''ADS test pattern.'''
    NORMAL = 0b000
    ALL_0 = 0b001
    ALL_1 = 0b010
    # 10101010101010 <=> 01010101010101
    TOGGLE = 0b011
    RAMP = 0b100
    CUSTOM = 0b101


class AdcSetting(RawSetting):
    '''Class to handle ADS4249 chip'''
    def __init__(self, rbcp_inst, verbose=True):
        RawSetting.__init__(self, rbcp_inst, verbose, 'ADC')
        self._flag_write = ADSWriteRead.read
        self.set_write_mode()

    def _write_reg(self, addr, data):
        self._write(ADS_OFFSET + addr, data)

    def _read_reg(self, addr):
        return self._read(ADS_OFFSET + addr)

    def reset(self):
        '''Reset the ADC chip'''
        self._write_reg(0x00, 0x02)
        self._flag_write = ADSWriteRead.write
        self.channel_swap(False)

    def _set_wr_mode(self, write_read):
        if self._flag_write == write_read:
            return

        self._flag_write = write_read
        self._write_reg(0x00, 0 if write_read == ADSWriteRead.write else 1)

    def set_write_mode(self):
        '''Configure the ADC chip for register writing'''
        self._set_wr_mode(ADSWriteRead.write)

    def set_read_mode(self):
        '''Configure the ADC chip to accept register reading'''
        self._set_wr_mode(ADSWriteRead.read)

    def write_reg(self, addr, data):
        '''Write to data to ADC registers

        Parameters
        ----------
        addr : int
            Register address.
        data : int
            Register value.
        '''
        self.set_write_mode()
        self._write_reg(addr, data)

    def read_reg(self, addr):
        '''Read ADC register.

        Parameters
        ----------
        addr : int
            Register address.

        Returns
        -------
        data : int
            Register data.
        '''
        self.set_read_mode()
        data = self._read_reg(addr)

        return data

    def _reg_ch_bit(self, addr, bit_num, val):
        # bit_num: 0-7
        # val: 0-1 or True/False
        val  = 1 if val else 0
        data = self.read_reg(addr)
        bit  = (data >> bit_num) & 0x01
        if bit == val:
            return

        mask1 = (1 << bit_num)
        mask2 = 0xff ^ mask1
        data = (mask1 * val) + (mask2 & data)
        self.write_reg(addr, data)

    def reg_ch_bit_up(self, addr, bit_num):
        '''Write 1 to the specified bit of the specified regsiter if it is currently 0.
        Do nothing when it is already 1.

        Parameters
        ----------
        addr : int
            Register address.
        bit_num : int
            Bit number.
        '''
        self._reg_ch_bit(addr, bit_num, True)


    def reg_ch_bit_down(self, addr, bit_num):
        '''Write 0 to the specified bit of the specified regsiter if it is currently 1.
        Do nothing when it is already 0.

        Parameters
        ----------
        addr : int
            Register address.
        bit_num : int
            Bit number.
        '''
        self._reg_ch_bit(addr, bit_num, False)

    @property
    def swapped(self):
        '''Tell whether channel definition is swapped.

        Parameters
        ----------
        swapped : bool
            True if swapped.
        '''
        return self._read(ADC_SWAP_ADDR) == 1

    def channel_swap(self, value=None):
        '''Swap channel A and B.

        Parameters
        ----------
        value : int, optional
            0...normal
            1...inverted
            If not specified, read the current value and then swap.
        '''
        if value is None:
            value = 0 if self.swapped  else 1 # swap
        else:
            value = 1 if value else 0

        assert isinstance(value, int)

        self._write(ADC_SWAP_ADDR, value)

    def set_digital(self, on_off):
        '''Turn on/off digital function.

        Parameters
        ----------
        on_off : bool
            True to turn on the digital function.
        '''

        if on_off:
            self.reg_ch_bit_up(0x42, 3)
        else:
            self.reg_ch_bit_down(0x42, 3)

    def set_digital_on(self):
        '''Turn on digital function.'''
        self.set_digital(True)

    def set_digital_off(self):
        '''Turn off digital function.'''
        self.set_digital(False)

    def set_test_pattern(self, channel:ADSChannel, pattern:ADSTestPattern):
        '''Set test pattern.

        Parameters
        ---------
        channel : ADSChannel, int
            A: 0, B: 1
        pattern : ADSTestPattern, int
            Test pattern.
        '''
        address = 0x25 if channel == ADSChannel.A else 0x2b

        pval = self.read_reg(address)
        cval = (pval & ~ADS_TEST_PATTERN_MASK) | pattern

        self.write_reg(address, cval)
