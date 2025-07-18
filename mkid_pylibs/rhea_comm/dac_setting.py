#!/usr/bin/env python3
# coding: utf-8
'''Handles DAC on the RHEA board'''

from raw_setting import RawSetting

DAC_OFFSET  = 0x20000000

DAC_TXEN    = 0x22000000
DAC_FRAME   = 0x22000001
DAC_SWAP    = 0x22000002
DAC_TESTEN  = 0x22000003
DAC_TESTPTN = 0x23000100

DAC_NORMAL_INPUT = 0
DAC_TEST_INPUT   = 1

ALARM_LABELS = ['alarm_fifo_1away',
                'alarm_fifo_2away',
                None,
                'alarm_from_iotest',
                None,
                'alarm_fifo_collision',
                'alarm_from_zerochk',
                None]

def parse_dac_alarm(reg_val):
    '''Parse register value of CONFIG7 of the DAC3283

    Parameter
    ---------
    reg_val : int
        CONFIG7 value.

    Returns
    -------
    alarm_dict : dict
        Dictionary that contains alarm information.
    '''
    alarm_dict = {}

    for i, label in enumerate(ALARM_LABELS):
        if label is None:
            continue

        mask = 1 << i
        alarm_dict[label] = bool(reg_val & mask)

    return alarm_dict


class DacSetting(RawSetting):
    '''Handles DAC3283 on the RHEA board.'''
    def __init__(self, rbcp_ins, verbose = True):
        RawSetting.__init__(self, rbcp_ins, verbose, 'DAC')
        self.flag_alarm = None
        self.flag_testmode = None
        self.en_4pin(force_reset = True)
        self.testmode_off(force_reset = True)

    def _write_reg(self, addr, data):
        self._write(DAC_OFFSET + addr, data)

    def _read_reg(self, addr):
        return self._read(DAC_OFFSET + addr)

    def en_4pin(self, force_reset=False):
        '''Enable 4-pin interface.

        Parameter
        ---------
        force_reset : bool, optional
            Force reset if True.
        '''
        if (not force_reset) and (not self.flag_alarm):
            return

        self.flag_alarm = False
        self._write_reg(0x17, 0x04)

    def en_alarm(self, force_reset=False):
        '''Enable alarm pin (instead of 4-pin interface).

        Parameter
        ---------
        force_reset : bool, optional
            Force reset if True.
        '''
        if (not force_reset) and self.flag_alarm:
            return

        self.flag_alarm = True
        self._write_reg(0x17, 0x00)

    def write_reg(self, addr, data):
        '''Write to a DAC register.

        Parameters
        ----------
        addr : int
            Register address.
        data : int
            Data to be wrote to the register.
        '''
        self._write_reg(addr, data)

    def read_reg(self, addr):
        '''Read register value.

        Parameter
        ---------
        addr : int
            Register address.
        '''
        self.en_4pin()
        return self._read_reg(addr)


    def reset(self):
        '''Reset DAC.
        The data sheet recommends to write default values to registers
        even though they are not being used.
        https://www.ti.com/jp/lit/ds/symlink/dac3283.pdf
        '''
        self.write_reg(0x00, 0x70)
        self.write_reg(0x01, 0x01) # Default 0x11 (FIR1 is on)
        self.write_reg(0x02, 0x00)
        self.write_reg(0x03, 0x13) # alarm_2/1away_ena
        self.write_reg(0x04, 0xff)
        self.write_reg(0x06, 0x00)
        self.write_reg(0x07, 0x00)
        self.write_reg(0x08, 0x00)
        self.write_reg(0x09, 0x00) # test pattern 0
        self.write_reg(0x0a, 0x00)
        self.write_reg(0x0b, 0x00)
        self.write_reg(0x0c, 0x00)
        self.write_reg(0x0d, 0x00)
        self.write_reg(0x0e, 0x00)
        self.write_reg(0x0f, 0x00)
        self.write_reg(0x10, 0x00) # test pattern 7
        self.write_reg(0x11, 0x24)
        self.write_reg(0x12, 0x02)
        self.write_reg(0x13, 0x02) # single_sync_source
        self.write_reg(0x14, 0x00)
        self.write_reg(0x15, 0x00)
        self.write_reg(0x16, 0x00)
        self.write_reg(0x17, 0x00)
        self.write_reg(0x14, 0x00) # again
        self.write_reg(0x18, 0x83)
        self.write_reg(0x19, 0x00)
        self.write_reg(0x1a, 0x00)
        self.write_reg(0x1b, 0x00)
        self.write_reg(0x1c, 0x00)
        self.write_reg(0x1d, 0x00)
        self.write_reg(0x1e, 0x24)
        self.write_reg(0x1b, 0x00) # again
        self.write_reg(0x1f, 0x12)
        self.en_4pin(force_reset=True)
        self.testmode_off(force_reset=True)

        self.txenable_off()
        self.trig_frame()
        self.channel_swap(False)
        self.dac_normal_input()

    def get_alarm(self):
        '''Get alarm information.

        Returns
        -------
        alarm_dict : dict
            Dictionary that contains alarm information.
        '''
        reg_val = self.read_reg(0x07)
        return parse_dac_alarm(reg_val)

    def clear_alarm(self):
        '''Clear alarm.'''
        self.write_reg(0x07, 0x00)

    def set_fifo_offset(self, offset):
        '''Change value of the fifo to control delay through the device.

        Parameter
        ---------
        offset : int
            Offset value [0, 8)
        '''
        assert 0 <= offset < 8
        cfg3 = self.read_reg(0x03)
        cfg3 &= 0xe3   # "11100011"
        cfg3 += (offset << 2)
        self.write_reg(0x03, cfg3)

    def get_fifo_offset(self):
        '''Read FIFO offset value.

        Returns
        -------
        offset : int
            FIFO offset value.
        '''
        return (self.read_reg(0x03) >> 2) & 0x07

    def set_testmode(self, on_off, force_reset=False):
        '''Turn on/off test mode.

        Parameters
        ----------
        on_off : bool
            Turn on if True.
        force_reset : bool, optional
            Force reset if True.
        '''
        if (not force_reset) and (self.flag_testmode == on_off):
            return

        self.flag_testmode = on_off
        cfg1 = self.read_reg(0x01)
        cfg1 = (cfg1 | 0x04) if on_off else (cfg1 & 0xfb)
        self.write_reg(0x01, cfg1)

    def testmode_on(self, force_reset=False):
        '''Trun on test mode.

        Parameter
        ---------
        force_reset : bool, optional
            Force reset if True.
        '''
        self.set_testmode(True, force_reset=force_reset)

    def testmode_off(self, force_reset = False):
        '''Trun off test mode.

        Parameter
        ---------
        force_reset : bool, optional
            Force reset if True.
        '''
        self.set_testmode(False, force_reset=force_reset)

    def set_test_pattern(self, adc_a, adc_b):
        '''Set test pattern to the DAC chip for communication test.

        Parameters
        ----------
        adc_a : int
            Value for Channel A.
            Should be less than 2^16
        adc_a : int
            Value for Channel B.
            Should be less than 2^16
        '''
        self.testmode_on()
        pt0 = (adc_a >> 8) & 0xff
        pt1 =  adc_a       & 0xff
        pt2 = (adc_b >> 8) & 0xff
        pt3 =  adc_b       & 0xff

        self.write_reg(0x09, pt0)
        self.write_reg(0x0a, pt1)
        self.write_reg(0x0b, pt2)
        self.write_reg(0x0c, pt3)
        self.write_reg(0x0d, pt0)
        self.write_reg(0x0e, pt1)
        self.write_reg(0x0f, pt2)
        self.write_reg(0x10, pt3)

    def clear_test_result(self):
        '''Clear test result.'''
        self.write_reg(0x08, 0x00)

    def get_test_result(self):
        '''Obtain iotest result.

        Returns
        -------
        iotest_results : int
            The values of these bits tell which bit in the byte-wide LVDS bus
            failed during the pattern checker test.
        '''

        return self.read_reg(0x08)

    def txenable_on(self):
        '''Set TXENABLE input to DAC high.'''
        self._write(DAC_TXEN, 1)

    def txenable_off(self):
        '''Set TXENABLE input to DAC low.'''
        self._write(DAC_TXEN, 0)

    def channel_swap(self, value=None):
        '''Swap Channel A and Channel B inside FPGA.

        Parameter
        ---------
        value : bool, optional
            Swap if True.
            If not specified, the function reads the current setting
            and change it to the opposite configuration.
        '''
        if value is None:
            data = self._read(DAC_SWAP)
            value = 0 if data else 1 # swap
        else:
            value = 1 if value else 0

        self._write(DAC_SWAP, value)

    def trig_frame(self):
        '''Trigger frame input.'''
        self._write(DAC_FRAME, 1)

    def dac_test_input(self, ch_a, ch_b, chmode=False):
        '''Set test input from FPGA to DAC for communication test.
        Use with `self.set_test_pattern` and `self.get_test_result`

        Parameters
        ----------
        ch_a : int
            Channel A test pattern.
        ch_b : int
            Channel B test pattern.
        chmode : bool, optional
            Change input source to test pattern if True.
        '''
        if chmode:
            self.set_input_source(DAC_TEST_INPUT)

        self._write2(DAC_TESTPTN, ch_a & 0xffff)
        self._write2(DAC_TESTPTN + 2, ch_b & 0xffff)

    def set_input_source(self, is_test):
        '''Set input source to the DAC chip.

        Parameter
        ---------
        is_test : int
            Set input source to test patterns if 1.
        '''
        self._write(DAC_TESTEN, is_test)

    def dac_normal_input(self):
        '''Set to the normal input mode.'''
        self.set_input_source(DAC_NORMAL_INPUT)
